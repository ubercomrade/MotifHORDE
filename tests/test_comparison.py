from __future__ import annotations

from dataclasses import replace

import pytest

from motifhorde.comparison import (
    compare_records,
    deduplicate_pairs,
    filter_similar,
    sort_pairs,
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


@pytest.mark.parametrize("metric", ["co", "dice", "cosine"])
def test_supported_profile_metrics_return_stable_pairs(pwm_model, sequences, metric):
    query = _record(pwm_model, "q")
    target = _record(replace(pwm_model, name="target"), "t")

    pairs = compare_records(query, [target], sequences, metric=metric, jobs=1)

    assert len(pairs) == 1
    assert pairs[0].query is query
    assert pairs[0].target is target
    assert pairs[0].result.metric == metric


def test_pair_filtering_sorting_and_deduplication(pwm_model, sequences):
    query = _record(pwm_model, "q")
    targets = [_record(replace(pwm_model, name=f"t{i}"), f"t{i}") for i in range(2)]
    pairs = compare_records(query, targets, sequences, metric="co", jobs=1)

    observed = deduplicate_pairs(
        sort_pairs(filter_similar(pairs, "score", -1), "score")
    )

    assert len(observed) == 1
    assert observed[0].target is targets[0]
