"""Mimosa profile comparison and stable list-based selection helpers."""

from __future__ import annotations

from dataclasses import dataclass, replace
from math import isfinite
from pathlib import Path
from typing import Any, Iterable

from mimosa import (
    EncodedSequences,
    HybridEmpiricalLogTail,
    MotifModel,
    NullDistribution,
    annotate_results,
    compare_many,
    normalization_fingerprint,
)
from mimosa.io import read_null_bundle, sequence_fingerprint

SUPPORTED_METRICS = ("co", "dice", "cosine")
DEFAULT_SCORE_COMPARISON_THRESHOLD = 0.9
DEFAULT_PVALUE_COMPARISON_THRESHOLD = 0.05
PROFILE_WINDOW_RADIUS = 10
PROFILE_REALIGN_WINDOW = 3
PROFILE_MIN_LOGERR = 0.0


@dataclass(frozen=True, slots=True)
class MatchedPair:
    """A Mimosa comparison result retaining the source records."""

    query: Any
    target: Any
    result: Any


@dataclass(frozen=True, slots=True)
class ComparisonSettings:
    """Validated Mimosa comparison parameters shared by all pipeline stages."""

    metric: str
    criterion: str
    threshold: float
    null_distribution: NullDistribution | None = None
    search_range: int = 10
    jobs: int = 1

    def __post_init__(self) -> None:
        if self.metric not in SUPPORTED_METRICS:
            raise ValueError(f"unsupported Mimosa metric: {self.metric}")
        if self.criterion not in {"score", "p-value"}:
            raise ValueError("comparison criterion must be 'score' or 'p-value'")
        if self.jobs < 1:
            raise ValueError("comparison jobs must be positive")
        if self.search_range < 0:
            raise ValueError("comparison search_range must be non-negative")
        if not isfinite(self.threshold):
            raise ValueError("comparison threshold must be finite")
        if self.criterion == "p-value" and self.null_distribution is None:
            raise ValueError("p-value comparison requires a Mimosa null bundle")

    def compare(
        self, query: Any, targets: Iterable[Any], sequences: EncodedSequences
    ) -> list[MatchedPair]:
        return compare_records(
            query,
            targets,
            sequences,
            metric=self.metric,
            jobs=self.jobs,
            null_distribution=self.null_distribution,
            search_range=self.search_range,
        )


def load_null_distribution(path: str | Path) -> NullDistribution:
    bundle = Path(path)
    if not bundle.is_dir():
        raise ValueError(f"Mimosa null distribution must be a directory bundle: {path}")
    return NullDistribution(**read_null_bundle(bundle))


def validate_null_distribution(
    distribution: NullDistribution,
    *,
    metric: str,
    search_range: int,
    sequences=None,
) -> None:
    """Check that a null bundle matches this module's public Mimosa call."""
    if not isinstance(distribution, NullDistribution):
        raise TypeError("null distribution must be a Mimosa NullDistribution")
    contract = distribution.contract
    expected_normalization = normalization_fingerprint(HybridEmpiricalLogTail())
    expected = {
        "metric": metric,
        "search_range": search_range,
        "window_radius": PROFILE_WINDOW_RADIUS,
        "realign_window": PROFILE_REALIGN_WINDOW,
        "min_logerr": PROFILE_MIN_LOGERR,
        "normalization_version": expected_normalization,
        "background_fingerprint": "none",
    }
    actual = {
        "metric": distribution.metric,
        "search_range": contract["search_range"],
        "window_radius": contract["window_radius"],
        "realign_window": contract["realign_window"],
        "min_logerr": float(contract["min_logerr"]),
        "normalization_version": contract["normalization_version"],
        "background_fingerprint": distribution.background_fingerprint,
    }
    if sequences is not None:
        expected["sequence_fingerprint"] = sequence_fingerprint(sequences)
        actual["sequence_fingerprint"] = distribution.sequence_fingerprint
    mismatches = [
        f"{key}: expected {expected[key]!r}, got {actual[key]!r}"
        for key in expected
        if actual[key] != expected[key]
    ]
    if mismatches:
        raise ValueError(
            "incompatible Mimosa null distribution: " + "; ".join(mismatches)
        )


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
    """Compare one query against targets in stable target order.

    The null distribution is used as given: callers passing a foreign
    bundle must validate it with validate_null_distribution themselves.
    """
    if metric not in SUPPORTED_METRICS:
        raise ValueError(f"unsupported Mimosa metric: {metric}")
    if jobs < 1:
        raise ValueError("comparison jobs must be positive")
    if search_range < 0:
        raise ValueError("comparison search_range must be non-negative")
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
        window_radius=PROFILE_WINDOW_RADIUS,
        realign_window=PROFILE_REALIGN_WINDOW,
        min_logerr=PROFILE_MIN_LOGERR,
        normalization=HybridEmpiricalLogTail(),
        total_threads=jobs,
        inner_threads=1,
    )
    if null_distribution is not None:
        results = annotate_results(
            results,
            null_distribution,
            effective_number_of_targets=len(target_records),
        )
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
