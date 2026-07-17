"""Motif scoring and binary classification metrics."""

function select_sequence_rows(batch::Mimosa.EncodedSequenceBatch, indices)
    rows = [copy(Mimosa.sequence(batch, Int(index))) for index in indices]
    return Mimosa.EncodedSequenceBatch(rows)
end

function all_valid_scores(model, sequences::Mimosa.EncodedSequenceBatch)
    scanned = Mimosa.scan(model, sequences; strands=Mimosa.BestStrand())
    return copy(scanned.data)
end

function best_scores(model, sequences::Mimosa.EncodedSequenceBatch)
    scanned = Mimosa.scan(model, sequences; strands=Mimosa.BestStrand())
    result = Float32[]
    for index in 1:Mimosa.nrows(scanned)
        values = Mimosa.row(scanned, index)
        isempty(values) || push!(result, maximum(values))
    end
    return result
end

struct PerformanceEvaluator
    background_type::Symbol
end
function PerformanceEvaluator(; background_type=:sites)
    return PerformanceEvaluator(Symbol(background_type))
end

function _roc_pr(labels::Vector{Bool}, scores::Vector{Float32})
    isempty(scores) && return (Float64[], Float64[], Float64[], Float64[], Float64[])
    order = sortperm(eachindex(scores); by=index -> (-scores[index], index))
    positives = count(labels)
    negatives = length(labels) - positives
    positives == 0 && return (Float64[], Float64[], Float64[], Float64[], Float64[])
    negatives == 0 && return (Float64[], Float64[], Float64[], Float64[], Float64[])
    recalls = Float64[0.0]
    precisions = Float64[1.0]
    fprs = Float64[0.0]
    tprs = Float64[0.0]
    thresholds = Float64[Inf]
    tp = 0
    fp = 0
    position = 1
    while position <= length(order)
        score = scores[order[position]]
        stop = position
        while stop <= length(order) && scores[order[stop]] == score
            labels[order[stop]] ? (tp += 1) : (fp += 1)
            stop += 1
        end
        push!(thresholds, Float64(score))
        push!(recalls, tp / positives)
        push!(precisions, tp + fp == 0 ? 1.0 : tp / (tp + fp))
        push!(fprs, fp / negatives)
        push!(tprs, tp / positives)
        position = stop
    end
    return recalls, precisions, fprs, tprs, thresholds
end

function _trapezoid(x, y)
    length(x) <= 1 && return 0.0
    return sum(
        (x[index + 1] - x[index]) * (y[index + 1] + y[index]) / 2 for
        index in 1:(length(x) - 1)
    )
end

function _standardized_pauc(area, minimum, maximum)
    denominator = maximum - minimum
    denominator == 0 && return 0.0
    return (area - minimum) / denominator
end

function evaluate(
    evaluator::PerformanceEvaluator, model, positives, negatives, error_threshold::Real
)
    positive_scores = best_scores(model, positives)
    negative_scores = if evaluator.background_type == :sites
        all_valid_scores(model, negatives)
    else
        best_scores(model, negatives)
    end
    labels = vcat(fill(true, length(positive_scores)), fill(false, length(negative_scores)))
    scores = vcat(positive_scores, negative_scores)
    recall, precision, fpr, tpr, thresholds = _roc_pr(labels, scores)
    auprc = _trapezoid(recall, precision)
    auroc = _trapezoid(fpr, tpr)
    table = Mimosa.fit(Mimosa.EmpiricalLogTail(), negative_scores)
    cutoff = Mimosa.lookup_score_for_tail_probability(table, Float64(error_threshold))
    roc_indices = findall(threshold -> threshold >= cutoff, thresholds)
    pr_indices = findall(threshold -> threshold >= cutoff, thresholds)
    cut_fpr = fpr[roc_indices]
    cut_tpr = tpr[roc_indices]
    cut_recall = recall[pr_indices]
    cut_precision = precision[pr_indices]
    endpoint_fpr = isempty(cut_fpr) ? 0.0 : cut_fpr[end]
    endpoint_recall = isempty(cut_recall) ? 0.0 : cut_recall[end]
    pauroc = _standardized_pauc(
        _trapezoid(cut_fpr, cut_tpr), endpoint_fpr^2 / 2, endpoint_fpr
    )
    pauprc = _standardized_pauc(
        _trapezoid(cut_recall, cut_precision), endpoint_recall / 2, endpoint_recall
    )
    return Dict{String,Any}(
        "PRC" => Dict("RECALL" => recall, "PRECISION" => precision),
        "ROC" => Dict("FPR" => fpr, "TPR" => tpr),
        "auPRC" => auprc,
        "auROC" => auroc,
        "pauPRC" => pauprc,
        "pauROC" => pauroc,
    )
end
