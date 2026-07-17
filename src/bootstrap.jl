"""Deterministic odd/even bootstrap task construction and execution."""

struct BootstrapTask
    index::Int
    params::Dict{String,Any}
    params_suffix::String
    split::Symbol
    foreground::String
    background::String
    output_dir::String
    number_of_motifs::Int
    seed::Union{Nothing,Int}
end

function _string_params(params)
    return Dict{String,Any}(String(key) => value for (key, value) in pairs(params))
end

function format_params(params::AbstractDict)
    normalized = _string_params(params)
    return join(
        ("$(key)-$(normalized[key])" for key in sort!(collect(keys(normalized)))), ","
    )
end

function parameter_grid(discovery_params)
    normalized = _string_params(discovery_params)
    key_names = sort!(collect(Base.keys(normalized)))
    values = [collect(normalized[key]) for key in key_names]
    return [
        Dict{String,Any}(key_names[i] => combination[i] for i in eachindex(key_names)) for
        combination in Iterators.product(values...)
    ]
end

function _bootstrap_indices(n::Int, split::Symbol)
    split in (:odd, :even) || throw(ArgumentError("Unknown bootstrap split: $split"))
    odd = [index for index in 1:n if isodd(index)]
    even = [index for index in 1:n if iseven(index)]
    return split == :odd ? (odd, even) : (even, odd)
end

function _task_seed(seed, index)
    return seed === nothing ? nothing : Int(seed) + index + 1
end

function build_bootstrap_tasks(
    peaks, background, discovery_params, task_root, number_of_motifs, seed
)
    tasks = BootstrapTask[]
    test_batches = Dict{Int,Mimosa.EncodedSequenceBatch}()
    index = 0
    for params in parameter_grid(discovery_params)
        suffix = format_params(params)
        for split in (:odd, :even)
            train_indices, test_indices = _bootstrap_indices(
                Mimosa.nsequences(peaks), split
            )
            task_dir = joinpath(
                task_root,
                "task_$(lpad(index, 4, '0'))_$(replace(suffix, '/' => '_'))_$(split)",
            )
            mkpath(task_dir)
            foreground_path = joinpath(task_dir, "train.fasta")
            background_path = joinpath(task_dir, "background.fasta")
            write_fasta(select_sequence_rows(peaks, train_indices), foreground_path)
            write_fasta(background, background_path)
            test_batches[index] = select_sequence_rows(peaks, test_indices)
            push!(
                tasks,
                BootstrapTask(
                    index,
                    params,
                    suffix,
                    split,
                    foreground_path,
                    background_path,
                    task_dir,
                    Int(number_of_motifs),
                    _task_seed(seed, index),
                ),
            )
            index += 1
        end
    end
    return tasks, test_batches
end

function _run_bootstrap_task(tool, task::BootstrapTask)
    params = copy(task.params)
    task.seed === nothing || (params["seed"] = task.seed)
    started = time()
    motifs = discover(
        tool,
        task.foreground,
        task.background,
        task.output_dir,
        task.number_of_motifs;
        _keyword_params(params)...,
    )
    return (
        index=task.index,
        params=task.params,
        params_suffix=task.params_suffix,
        split=task.split,
        elapsed=time() - started,
        motifs=motifs,
    )
end

function _run_tasks(tool, tasks, jobs)
    results = Vector{Any}(undef, length(tasks))
    if jobs <= 1 || Threads.nthreads() == 1
        for index in eachindex(tasks)
            results[index] = _run_bootstrap_task(tool, tasks[index])
        end
    else
        Threads.@threads for index in eachindex(tasks)
            results[index] = _run_bootstrap_task(tool, tasks[index])
        end
    end
    return results
end

function run_bootstrap(
    tool,
    evaluator,
    output_dir,
    peaks,
    background,
    number_of_motifs,
    error_threshold,
    discovery_params;
    jobs=1,
    seed=nothing,
)
    task_parent = joinpath(output_dir, "bootstrap")
    mkpath(task_parent)
    task_root = mktempdir(task_parent; prefix="bootstrap_")
    tasks, test_batches = build_bootstrap_tasks(
        peaks, background, discovery_params, task_root, number_of_motifs, seed
    )
    results = sort!(_run_tasks(tool, tasks, max(1, Int(jobs))); by=result -> result.index)
    statistics = Dict{String,Any}()
    motifs = Mimosa.AbstractMotifModel[]
    for result in results
        for model in result.motifs
            model_stats = evaluate(
                evaluator, model, test_batches[result.index], background, error_threshold
            )
            name = "$(motif_name(model))_$(result.params_suffix)_$(result.split)"
            renamed = rename_model(model, name)
            statistics[name] = model_stats
            push!(motifs, renamed)
        end
    end
    return statistics, motifs
end
