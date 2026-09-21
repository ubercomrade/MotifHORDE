"""Mimosa profile comparison and stable list-based selection helpers."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Iterable

from mimosa import (
    MotifModel,
    NullDistribution,
    annotate_results,
    compare_many,
)
from mimosa.io import read_null_bundle

SUPPORTED_METRICS = ("co", "dice", "cosine")
DEFAULT_SCORE_COMPARISON_THRESHOLD = 0.9
DEFAULT_PVALUE_COMPARISON_THRESHOLD = 0.05


@dataclass(frozen=True, slots=True)
class MatchedPair:
    """A Mimosa comparison result retaining the source records."""

    query: Any
    target: Any
    result: Any


def load_null_distribution(path: str | Path) -> NullDistribution:
    bundle = Path(path)
    if not bundle.is_dir():
        raise ValueError(f"Mimosa null distribution must be a directory bundle: {path}")
    return read_null_bundle(bundle)


def _model(record: Any) -> MotifModel:
    model = getattr(record, "model", record)
    if not isinstance(model, MotifModel):
        raise TypeError(f"comparison input is not a MotifModel: {type(model).__name__}")
    return model


def _identity(record: Any, model: MotifModel) -> str:
    return str(getattr(record, "id", model.name))


def _named(model: MotifModel, name: str) -> MotifModel:
    return replace(model, name=name)


def compare_records(
    query: Any,
    targets: Iterable[Any],
    sequences,
    *,
    metric: str,
    jobs: int,
    null_distribution: NullDistribution | None = None,
    search_range: int = 10,
) -> list[MatchedPair]:
    """Compare one query against targets in stable target order."""
    if metric not in SUPPORTED_METRICS:
        raise ValueError(f"unsupported Mimosa metric: {metric}")
    if jobs < 1:
        raise ValueError("comparison jobs must be positive")
    target_records = list(targets)
    query_model = _model(query)
    query_for_mimosa = _named(query_model, _identity(query, query_model))
    target_models = [
        _named(target_model, _identity(record, target_model))
        for record in target_records
        for target_model in [_model(record)]
    ]
    results = compare_many(
        query_for_mimosa,
        target_models,
        sequences,
        metric=metric,
        search_range=search_range,
        total_threads=jobs,
        inner_threads=1,
    )
    if null_distribution is not None:
        results = annotate_results(results, null_distribution)
    return [
        MatchedPair(query, target, result)
        for target, result in zip(target_records, results, strict=True)
    ]


def comparison_value(pair: MatchedPair, criterion: str) -> float:
    if criterion == "score":
        return float(pair.result.score)
    if criterion == "p-value":
        value = getattr(pair.result, "adj_p_value", None)
        if value is None:
            raise ValueError(
                "Mimosa p-value comparison requires a compatible null distribution bundle"
            )
        return float(value)
    raise ValueError(f"unsupported comparison criterion: {criterion}")


def default_threshold(criterion: str) -> float:
    if criterion == "score":
        return DEFAULT_SCORE_COMPARISON_THRESHOLD
    if criterion == "p-value":
        return DEFAULT_PVALUE_COMPARISON_THRESHOLD
    raise ValueError(f"unsupported comparison criterion: {criterion}")


def is_similar(pair: MatchedPair, criterion: str, threshold: float) -> bool:
    value = comparison_value(pair, criterion)
    return value >= threshold if criterion == "score" else value <= threshold


def filter_similar(
    pairs: Iterable[MatchedPair], criterion: str, threshold: float
) -> list[MatchedPair]:
    return [pair for pair in pairs if is_similar(pair, criterion, threshold)]


def sort_pairs(pairs: Iterable[MatchedPair], criterion: str) -> list[MatchedPair]:
    return sorted(
        pairs,
        key=lambda pair: comparison_value(pair, criterion),
        reverse=criterion == "score",
    )


def deduplicate_pairs(pairs: Iterable[MatchedPair]) -> list[MatchedPair]:
    """Keep one stable query/target matching from already sorted pairs."""
    result = []
    queries: set[str] = set()
    targets: set[str] = set()
    for pair in pairs:
        query_model = _model(pair.query)
        target_model = _model(pair.target)
        query_id = _identity(pair.query, query_model)
        target_id = _identity(pair.target, target_model)
        if query_id in queries or target_id in targets:
            continue
        queries.add(query_id)
        targets.add(target_id)
        result.append(pair)
    return result
