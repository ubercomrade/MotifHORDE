from __future__ import annotations

import pandas as pd
import pytest

from motifhorde.models import GenericModel
from motifhorde.pipeline import (
    DeNovoPipeline,
    _comparison_column,
    _deduplicate_final_motifs,
    _deduplicate_matches,
    _filter_similar_matches,
    _is_similar_value,
    _select_nonredundant_motifs,
    _sort_comparisons,
)


def motif(name: str) -> GenericModel:
    return GenericModel("test", name, None, 1, {})


class FakeComparator:
    def __init__(self, column: str, values: dict[tuple[str, str], float]) -> None:
        self.column = column
        self.values = values

    def compare(self, motifs_1, motifs_2, sequences=None):
        return pd.DataFrame(
            [
                {"query": query.name, "target": target.name, self.column: value}
                for query in motifs_1
                for target in motifs_2
                if (value := self.values.get((query.name, target.name))) is not None
            ],
            columns=["query", "target", self.column],
        )


def stats(names: list[str], metric_values: list[float]) -> dict[str, dict[str, float]]:
    return {
        name: {
            "auPRC": value,
            "auROC": value,
            "pauPRC": value,
            "pauROC": value,
        }
        for name, value in zip(names, metric_values)
    }


def test_similarity_thresholds_are_inclusive():
    assert _is_similar_value("p-value", 0.001)
    assert not _is_similar_value("p-value", 0.0011)
    assert _is_similar_value("score", 0.9)
    assert not _is_similar_value("score", 0.899)


def test_sort_comparisons_uses_comparison_direction():
    pvalue_frame = pd.DataFrame(
        [
            {"query": "q2", "target": "t2", "p-value": 0.002},
            {"query": "q1", "target": "t1", "p-value": 0.001},
        ]
    )
    score_frame = pd.DataFrame(
        [
            {"query": "q2", "target": "t2", "score": 0.91},
            {"query": "q1", "target": "t1", "score": 0.95},
        ]
    )

    assert list(_sort_comparisons(pvalue_frame)["query"]) == ["q1", "q2"]
    assert list(_sort_comparisons(score_frame)["query"]) == ["q1", "q2"]

    with pytest.raises(ValueError):
        _comparison_column(pd.DataFrame([{"query": "q", "target": "t"}]))


def test_deduplicate_matches_keeps_first_query_and_target_matches():
    frame = pd.DataFrame(
        [
            {"query": "q1", "target": "t1", "p-value": 0.0001},
            {"query": "q1", "target": "t2", "p-value": 0.0002},
            {"query": "q2", "target": "t1", "p-value": 0.0003},
            {"query": "q3", "target": "t3", "p-value": 0.0004},
        ]
    )

    observed = _deduplicate_matches(frame)

    assert list(observed["query"]) == ["q1", "q3"]
    assert list(observed["target"]) == ["t1", "t3"]


def test_select_nonredundant_motifs_keeps_best_metric_representative():
    motifs = [motif("low"), motif("high"), motif("distinct")]
    comparator = FakeComparator(
        "score",
        {
            ("low", "high"): 0.95,
            ("distinct", "high"): 0.2,
        },
    )

    selected = _select_nonredundant_motifs(
        motifs,
        stats(["low", "high", "distinct"], [0.5, 0.9, 0.4]),
        "pauROC",
        comparator,
        sequences=None,
    )

    assert [model.name for model in selected] == ["high", "distinct"]


def test_bootstrap_matching_is_parameter_local():
    bootstrap_motifs = [
        motif("m1_length-8_odd"),
        motif("m2_length-8_odd"),
        motif("m3_length-8_even"),
        motif("m4_length-10_odd"),
        motif("m5_length-10_even"),
    ]
    comparator = FakeComparator(
        "score",
        {
            ("m1_length-8_odd", "m2_length-8_odd"): 0.95,
            ("m1_length-8_odd", "m3_length-8_even"): 0.96,
            ("m4_length-10_odd", "m5_length-10_even"): 0.97,
            ("m1_length-8_odd", "m5_length-10_even"): 0.99,
        },
    )
    pipeline = DeNovoPipeline(None, None, comparator)

    frame = pipeline._compare_bootstrap_motifs(
        bootstrap_motifs,
        {"length": [8, 10]},
        peaks=None,
        statistics=stats([model.name for model in bootstrap_motifs], [0.9, 0.8, 0.7, 0.6, 0.5]),
        metric="pauROC",
    )

    assert frame is not None
    assert set(zip(frame["query"], frame["target"])) == {
        ("m1_length-8_odd", "m3_length-8_even"),
        ("m4_length-10_odd", "m5_length-10_even"),
    }
    assert set(frame["length"]) == {8, 10}


def test_select_best_full_motif_requires_similarity_to_both_references_pvalue():
    full_motifs = [motif("full_odd_only"), motif("full_even_only"), motif("full_best"), motif("full_second")]
    bootstrap_motifs = [motif("odd_ref"), motif("even_ref")]
    comparator = FakeComparator(
        "p-value",
        {
            ("full_odd_only", "odd_ref"): 0.0001,
            ("full_odd_only", "even_ref"): 0.01,
            ("full_even_only", "odd_ref"): 0.01,
            ("full_even_only", "even_ref"): 0.0001,
            ("full_best", "odd_ref"): 0.0002,
            ("full_best", "even_ref"): 0.0004,
            ("full_second", "odd_ref"): 0.0005,
            ("full_second", "even_ref"): 0.0005,
        },
    )
    pipeline = DeNovoPipeline(None, None, comparator)

    selected = pipeline._select_best_full_motif(full_motifs, "odd_ref", "even_ref", bootstrap_motifs, peaks=None)

    assert selected is not None
    assert selected.name == "full_best"


def test_select_best_full_motif_uses_highest_average_score():
    full_motifs = [motif("full_low"), motif("full_best")]
    bootstrap_motifs = [motif("odd_ref"), motif("even_ref")]
    comparator = FakeComparator(
        "score",
        {
            ("full_low", "odd_ref"): 0.99,
            ("full_low", "even_ref"): 0.90,
            ("full_best", "odd_ref"): 0.96,
            ("full_best", "even_ref"): 0.96,
        },
    )
    pipeline = DeNovoPipeline(None, None, comparator)

    selected = pipeline._select_best_full_motif(full_motifs, "odd_ref", "even_ref", bootstrap_motifs, peaks=None)

    assert selected is not None
    assert selected.name == "full_best"


def test_filter_similar_matches_rejects_non_similar_rows():
    frame = pd.DataFrame(
        [
            {"query": "q1", "target": "t1", "score": 0.9},
            {"query": "q2", "target": "t2", "score": 0.899},
        ]
    )

    observed = _filter_similar_matches(frame)

    assert list(observed["query"]) == ["q1"]


def test_final_global_deduplication_keeps_best_metric_motif():
    final_motifs = [motif("weak"), motif("strong"), motif("distinct")]
    final_info = [
        ("weak", {"length": 8}),
        ("strong", {"length": 10}),
        ("distinct", {"length": 12}),
    ]
    final_stats = {
        "weak_length-8": {"pauROC": 0.7},
        "strong_length-10": {"pauROC": 0.9},
        "distinct_length-12": {"pauROC": 0.8},
    }
    comparator = FakeComparator(
        "score",
        {
            ("weak", "strong"): 0.95,
            ("distinct", "strong"): 0.1,
        },
    )

    motifs, info, stats_observed = _deduplicate_final_motifs(
        final_motifs,
        final_info,
        final_stats,
        "pauROC",
        comparator,
        sequences=None,
    )

    assert [model.name for model in motifs] == ["strong", "distinct"]
    assert info == [("strong", {"length": 10}), ("distinct", {"length": 12})]
    assert set(stats_observed) == {"strong_length-10", "distinct_length-12"}
