from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from mimosa import build_null
from mimosa.io import write_null_bundle

from motifhorde.comparison import (
    compare_records,
    comparison_sequences,
    deduplicate_pairs,
    filter_similar,
    is_similar,
    load_null_distribution,
    sort_pairs,
    validate_null_distribution,
)
from motifhorde.evaluation import EvaluatedMotif


def _record(model, identifier, metrics=None):
    return EvaluatedMotif(
        model,
        identifier,
        identifier.split(":")[0],
        "odd",
        {},
        metrics or {"auPRC": 0.5, "auROC": 0.5, "pauPRC": 0.5, "pauROC": 0.5},
    )


def test_profile_comparison_returns_stable_jaccard_pairs(pwm_model, sequences):
    query = _record(pwm_model, "q")
    target = _record(replace(pwm_model, name="target"), "t")

    pairs = compare_records(query, [target], sequences, jobs=1)

    assert len(pairs) == 1
    assert pairs[0].query is query
    assert pairs[0].target is target
    assert pairs[0].result.metric == "jaccard"


def test_comparison_sequences_are_fixed_random_set():
    first = comparison_sequences()
    second = comparison_sequences()

    assert len(first) == 10_000
    assert set(np.diff(first.offsets).tolist()) == {100}
    assert set(np.unique(first.data).tolist()) == {0, 1, 2, 3}
    assert first == second


def test_pair_filtering_sorting_and_deduplication(pwm_model, sequences):
    query = _record(pwm_model, "q")
    targets = [_record(replace(pwm_model, name=f"t{i}"), f"t{i}") for i in range(2)]
    pairs = compare_records(query, targets, sequences, jobs=1)

    observed = deduplicate_pairs(
        sort_pairs(filter_similar(pairs, "score", -1), "score")
    )

    assert len(observed) == 1
    assert observed[0].target is targets[0]


def _write_null_bundle(tmp_path: Path, pwm_model, sequences, *, search_range=10):
    other = replace(pwm_model, name="Other")
    distribution = build_null(
        [pwm_model, other],
        sequences=sequences,
        n_samples=16,
        seed=7,
        search_range=search_range,
    )
    path = tmp_path / "null"
    write_null_bundle(path, distribution)
    return path


def test_pvalue_comparison_loads_directory_bundle_and_applies_bh_per_query(
    tmp_path, pwm_model, sequences, monkeypatch
):
    path = _write_null_bundle(tmp_path, pwm_model, sequences)
    null_distribution = load_null_distribution(path)
    query = _record(pwm_model, "query")
    targets = [
        _record(replace(pwm_model, name=f"target-{index}"), f"target-{index}")
        for index in range(2)
    ]
    calls = []

    import motifhorde.comparison as comparison

    original = comparison.annotate_results

    def annotate(results, distribution, **kwargs):
        calls.append(kwargs["effective_number_of_targets"])
        return original(results, distribution, **kwargs)

    monkeypatch.setattr(comparison, "annotate_results", annotate)
    pairs = compare_records(
        query,
        targets,
        sequences,
        jobs=1,
        null_distribution=null_distribution,
    )
    compare_records(
        _record(replace(pwm_model, name="query-2"), "query-2"),
        [targets[0]],
        sequences,
        jobs=1,
        null_distribution=null_distribution,
    )

    assert calls == [2, 1]
    assert all(0 <= pair.result.p_value <= 1 for pair in pairs)
    assert all(0 <= pair.result.adj_p_value <= 1 for pair in pairs)


def test_null_bundle_validation_rejects_incompatible_search_range(
    tmp_path, pwm_model, sequences
):
    path = _write_null_bundle(tmp_path, pwm_model, sequences, search_range=9)
    null_distribution = load_null_distribution(path)

    with pytest.raises(ValueError, match="incompatible"):
        validate_null_distribution(
            null_distribution, search_range=10, sequences=sequences
        )
    validate_null_distribution(null_distribution, search_range=9, sequences=sequences)


def test_null_bundle_validation_rejects_other_min_logerr(
    tmp_path, pwm_model, sequences
):
    distribution = build_null(
        [pwm_model, replace(pwm_model, name="Other")],
        sequences=sequences,
        n_samples=16,
        seed=7,
        search_range=10,
        min_logerr=2.0,
    )
    path = tmp_path / "null"
    write_null_bundle(path, distribution)

    with pytest.raises(ValueError, match="min_logerr"):
        validate_null_distribution(
            load_null_distribution(path), search_range=10, sequences=sequences
        )


def test_corrupt_null_bundle_is_rejected(tmp_path):
    bundle = tmp_path / "null"
    bundle.mkdir()
    (bundle / "manifest.toml").write_text("not a Mimosa bundle")

    with pytest.raises(Exception):
        load_null_distribution(bundle)


def test_thresholds_are_inclusive():
    pair = SimpleNamespace(result=SimpleNamespace(score=0.9))
    assert is_similar(pair, "score", 0.9)
    pair.result = SimpleNamespace(adj_p_value=0.05)
    assert is_similar(pair, "p-value", 0.05)
