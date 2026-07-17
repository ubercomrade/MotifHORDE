"""Pipeline orchestration; computation remains in small testable functions."""

const VALIDATION_METRICS = ("auPRC", "auROC", "pauPRC", "pauROC")

struct PipelineConfig
    fpr_threshold::Float64
    number_of_motifs::Int
    jobs::Int
    seed::Union{Nothing,Int}
end

function _tool_name(tool::MotifDiscoveryTool)
    names = Dict(
        StremeDiscoveryTool => "streme",
        MemeDiscoveryTool => "meme",
        BammDiscoveryTool => "bamm",
        DimontDiscoveryTool => "dimont",
        SlimDiscoveryTool => "slim",
        SitegaDiscoveryTool => "sitega",
    )
    return get(names, typeof(tool), lowercase(string(nameof(typeof(tool)))))
end

function _keyword_params(params)
    names = Tuple(Symbol.(collect(keys(params))))
    return NamedTuple{names}(Tuple(params[String(name)] for name in names))
end

function _bootstrap_motifs_for_params(models, params)
    suffix = format_params(params)
    odd_suffix = "_$(suffix)_odd"
    even_suffix = "_$(suffix)_even"
    odd = [model for model in models if endswith(motif_name(model), odd_suffix)]
    even = [model for model in models if endswith(motif_name(model), even_suffix)]
    return odd, even
end

function _compare_bootstrap(models, params_grid, statistics, comparator, sequences, metric)
    records = NamedTuple[]
    for params in params_grid
        odd, even = _bootstrap_motifs_for_params(models, params)
        (isempty(odd) || isempty(even)) && continue
        selected_odd = select_nonredundant_motifs(
            odd, statistics, metric, comparator, sequences
        )
        selected_even = select_nonredundant_motifs(
            even, statistics, metric, comparator, sequences
        )
        rows = Comparison[]
        for model in selected_odd
            append!(
                rows,
                filter_similar(
                    compare(model, selected_even, sequences, comparator), comparator
                ),
            )
        end
        for row in deduplicate_matches(sort_comparisons(rows, comparator))
            query_stats = get(statistics, row.query, Dict{String,Any}())
            target_stats = get(statistics, row.target, Dict{String,Any}())
            validation = Dict{String,Float64}(
                validation_metric =>
                    (
                        Float64(get(query_stats, validation_metric, 0.0)) +
                        Float64(get(target_stats, validation_metric, 0.0))
                    ) / 2 for validation_metric in VALIDATION_METRICS
            )
            push!(
                records,
                (
                    params=params,
                    query=row.query,
                    target=row.target,
                    score=row.score,
                    adjusted_p_value=row.adjusted_p_value,
                    validation=validation,
                ),
            )
        end
    end
    return records
end

function _best_full_model(
    full_models, odd_name, even_name, bootstrap_models, comparator, sequences
)
    odd = findfirst(model -> motif_name(model) == odd_name, bootstrap_models)
    even = findfirst(model -> motif_name(model) == even_name, bootstrap_models)
    (odd === nothing || even === nothing) && return nothing
    candidates = Tuple{Any,Float64}[]
    for model in full_models
        odd_row = only_or_nothing(
            compare(model, [bootstrap_models[odd]], sequences, comparator)
        )
        even_row = only_or_nothing(
            compare(model, [bootstrap_models[even]], sequences, comparator)
        )
        if odd_row === nothing || even_row === nothing
            continue
        end
        _is_similar_value(comparator, _comparison_value(odd_row, comparator)) || continue
        _is_similar_value(comparator, _comparison_value(even_row, comparator)) || continue
        average =
            (
                _comparison_value(odd_row, comparator) +
                _comparison_value(even_row, comparator)
            ) / 2
        push!(candidates, (model, average))
    end
    isempty(candidates) && return nothing
    return sort(
        candidates; by=pair -> pair[2], rev=comparator.comparison_criterion == :score
    )[1][1]
end

only_or_nothing(values) = isempty(values) ? nothing : first(values)

function _ensure_pfm(model, promoters)
    return Mimosa.reconstruct_pfm(
        model, promoters; mode=:best, strands=Mimosa.BothStrands()
    )
end

function _write_json(path, value)
    open(path, "w") do io
        return JSON3.write(io, value; indent=2)
    end
    return path
end

function _save_models(models, statistics, directory, promoters)
    models_dir = joinpath(directory, "models")
    mkpath(models_dir)
    pfms = Matrix{Float32}[]
    metadata = Tuple{String,Int}[]
    for (rank, model) in enumerate(models)
        model_path = joinpath(
            models_dir, @sprintf("%03d_%s", rank, replace(motif_name(model), '/' => '_'))
        )
        write_model(model_path, model)
        pfm = _ensure_pfm(model, promoters)
        push!(pfms, pfm)
        push!(metadata, (motif_name(model), motif_length(model)))
    end
    write_meme(pfms, metadata, joinpath(models_dir, "all_motifs_in_pfm_form.meme"))
    _write_json(joinpath(directory, "statistics.json"), statistics)
    return directory
end

function run_pipeline(
    tool::MotifDiscoveryTool,
    evaluator::PerformanceEvaluator,
    comparator::AbstractComparator,
    foreground_path::AbstractString,
    background_path::AbstractString,
    promoters_path::AbstractString,
    output_dir::AbstractString,
    discovery_params;
    metric::AbstractString="pauROC",
    config::PipelineConfig=PipelineConfig(0.001, 5, 1, nothing),
)
    metric in VALIDATION_METRICS ||
        throw(ArgumentError("Unsupported validation metric: $metric"))
    all(isfile, (foreground_path, background_path, promoters_path)) ||
        throw(ArgumentError("All FASTA input paths must exist."))
    foreground, _ = read_fasta(foreground_path)
    background, _ = read_fasta(background_path)
    promoters, _ = read_fasta(promoters_path)
    tool_dir = joinpath(output_dir, _tool_name(tool))
    mkpath(tool_dir)
    statistics, bootstrap_models = run_bootstrap(
        tool,
        evaluator,
        tool_dir,
        foreground,
        background,
        config.number_of_motifs,
        config.fpr_threshold,
        discovery_params;
        jobs=config.jobs,
        seed=config.seed,
    )
    bootstrap_dir = joinpath(tool_dir, "bootstrap")
    _save_models(bootstrap_models, statistics, bootstrap_dir, promoters)
    params_grid = parameter_grid(discovery_params)
    records = _compare_bootstrap(
        bootstrap_models, params_grid, statistics, comparator, foreground, metric
    )
    final_models = Mimosa.AbstractMotifModel[]
    final_info = Tuple{String,Dict{String,Any}}[]
    final_statistics = Dict{String,Any}()
    for record in records
        params = record.params
        full_dir = mktempdir(tool_dir; prefix="final_")
        full_models = discover(
            tool,
            foreground_path,
            background_path,
            full_dir,
            config.number_of_motifs * 2;
            _keyword_params(params)...,
        )
        best = _best_full_model(
            full_models,
            record.query,
            record.target,
            bootstrap_models,
            comparator,
            foreground,
        )
        best === nothing && continue
        key = "$(motif_name(best))_$(format_params(params))"
        push!(final_models, best)
        push!(final_info, (motif_name(best), params))
        final_statistics[key] = record.validation
    end
    final_models, final_info, final_statistics = deduplicate_final_motifs(
        final_models, final_info, final_statistics, metric, comparator, foreground
    )
    final_dir = joinpath(tool_dir, "motifs")
    _save_models(final_models, final_statistics, final_dir, promoters)
    return final_models
end
