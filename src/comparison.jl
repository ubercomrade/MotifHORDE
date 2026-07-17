"""Mimosa-backed comparison and pure selection rules."""

const DEFAULT_SCORE_COMPARISON_THRESHOLD = 0.9
const DEFAULT_PVALUE_COMPARISON_THRESHOLD = 0.05
const MISSING_ADJUSTED_PVALUE_MESSAGE = "Comparison results do not contain adjusted p-values. Provide a compatible Mimosa null distribution."

struct Comparison
    query::String
    target::String
    score::Float64
    adjusted_p_value::Union{Nothing,Float64}
    offset::Int
    orientation::String
end

abstract type AbstractComparator end

struct TomtomComparator <: AbstractComparator
    metric::Symbol
    comparison_criterion::Symbol
    comparison_threshold::Float64
    null_distribution::Union{Nothing,String}
end

function TomtomComparator(;
    metric=:pcc,
    comparison_criterion=:score,
    comparison_threshold=nothing,
    null_distribution=nothing,
    kwargs...,
)
    criterion = Symbol(comparison_criterion)
    threshold = if comparison_threshold === nothing
        default_threshold_for_criterion(criterion)
    else
        Float64(comparison_threshold)
    end
    return TomtomComparator(
        Symbol(metric),
        criterion,
        threshold,
        null_distribution === nothing ? nothing : String(null_distribution),
    )
end

struct UniversalMotifComparator <: AbstractComparator
    metric::Symbol
    comparison_criterion::Symbol
    comparison_threshold::Float64
    search_range::Int
    null_distribution::Union{Nothing,String}
end

function UniversalMotifComparator(;
    metric=:co,
    comparison_criterion=:score,
    comparison_threshold=nothing,
    search_range=10,
    null_distribution=nothing,
    kwargs...,
)
    criterion = Symbol(comparison_criterion)
    threshold = if comparison_threshold === nothing
        default_threshold_for_criterion(criterion)
    else
        Float64(comparison_threshold)
    end
    return UniversalMotifComparator(
        Symbol(metric),
        criterion,
        threshold,
        Int(search_range),
        null_distribution === nothing ? nothing : String(null_distribution),
    )
end

function comparison_column_for_criterion(criterion::Symbol)
    if criterion == :score
        :score
    elseif criterion == Symbol("p-value")
        Symbol("adj.p-value")
    else
        throw(ArgumentError("Unsupported comparison criterion: $criterion"))
    end
end
function comparison_column_for_criterion(criterion::AbstractString)
    return comparison_column_for_criterion(Symbol(criterion))
end

function default_threshold_for_criterion(criterion::Symbol)
    if criterion == :score
        DEFAULT_SCORE_COMPARISON_THRESHOLD
    elseif criterion == Symbol("p-value")
        DEFAULT_PVALUE_COMPARISON_THRESHOLD
    else
        throw(ArgumentError("Unsupported comparison criterion: $criterion"))
    end
end
function default_threshold_for_criterion(criterion::AbstractString)
    return default_threshold_for_criterion(Symbol(criterion))
end

function _raw_comparison(query, target, sequences, comparator)
    sequences === nothing &&
        throw(ArgumentError("Mimosa profile comparison requires sequences."))
    metric = comparator isa UniversalMotifComparator ? comparator.metric : :co
    search_range = comparator isa UniversalMotifComparator ? comparator.search_range : 10
    return Mimosa.compare(
        query, target, sequences; metric=metric, search_range=search_range
    )
end

function _comparison_rows(query, targets, sequences, comparator)
    raw = Mimosa.ComparisonResult[]
    for target in targets
        push!(raw, _raw_comparison(query, target, sequences, comparator))
    end
    null_distribution = if comparator isa TomtomComparator
        comparator.null_distribution
    else
        comparator.null_distribution
    end
    if null_distribution !== nothing
        distribution = Mimosa.loadnull(null_distribution)
        annotated = Mimosa.annotate_results(raw, distribution)
        return [
            Comparison(
                r.query, r.target, Float64(r.score), r.adj_p_value, r.offset, r.orientation
            ) for r in annotated
        ]
    end
    return [
        Comparison(r.query, r.target, Float64(r.score), nothing, r.offset, r.orientation)
        for r in raw
    ]
end

function compare(query, targets, sequences, comparator::AbstractComparator)
    return _comparison_rows(query, targets, sequences, comparator)
end

function _comparison_value(row::Comparison, comparator)
    comparator.comparison_criterion == :score && return row.score
    row.adjusted_p_value === nothing &&
        throw(ArgumentError(MISSING_ADJUSTED_PVALUE_MESSAGE))
    return row.adjusted_p_value::Float64
end

function _is_similar_value(comparator::AbstractComparator, value)
    isfinite(Float64(value)) || return false
    comparator.comparison_criterion == :score &&
        return value >= comparator.comparison_threshold
    comparator.comparison_criterion == Symbol("p-value") &&
        return value <= comparator.comparison_threshold
    return throw(
        ArgumentError(
            "Unsupported comparison criterion: $(comparator.comparison_criterion)"
        ),
    )
end

function filter_similar(rows::AbstractVector{Comparison}, comparator::AbstractComparator)
    return [
        row for
        row in rows if _is_similar_value(comparator, _comparison_value(row, comparator))
    ]
end

function sort_comparisons(rows::AbstractVector{Comparison}, comparator::AbstractComparator)
    return sort(
        rows;
        by=row -> _comparison_value(row, comparator),
        rev=comparator.comparison_criterion == :score,
    )
end

function deduplicate_matches(rows::AbstractVector{Comparison})
    seen_query = Set{String}()
    seen_target = Set{String}()
    result = Comparison[]
    for row in rows
        row.query in seen_query && continue
        row.target in seen_target && continue
        push!(result, row)
        push!(seen_query, row.query)
        push!(seen_target, row.target)
    end
    return result
end

function select_nonredundant_motifs(models, statistics, metric, comparator, sequences)
    ordered = sort(
        collect(models);
        by=model -> get(get(statistics, motif_name(model), Dict()), metric, 0.0),
        rev=true,
    )
    selected = Mimosa.AbstractMotifModel[]
    for model in ordered
        isempty(selected) && (push!(selected, model); continue)
        rows = filter_similar(compare(model, selected, sequences, comparator), comparator)
        isempty(rows) && push!(selected, model)
    end
    return selected
end

function deduplicate_final_motifs(models, info, statistics, metric, comparator, sequences)
    order = sortperm(
        eachindex(models);
        by=index -> get(
            statistics,
            "$(info[index][1])_$(format_params(info[index][2]))",
            Dict(metric => 0.0),
        )[metric],
        rev=true,
    )
    kept = Mimosa.AbstractMotifModel[]
    kept_info = Tuple{String,Dict{String,Any}}[]
    kept_stats = Dict{String,Dict{String,Float64}}()
    for index in order
        rows = if isempty(kept)
            Comparison[]
        else
            filter_similar(compare(models[index], kept, sequences, comparator), comparator)
        end
        isempty(rows) || continue
        push!(kept, models[index])
        push!(kept_info, info[index])
        key = "$(info[index][1])_$(format_params(info[index][2]))"
        kept_stats[key] = statistics[key]
    end
    return kept, kept_info, kept_stats
end
