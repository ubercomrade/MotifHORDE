from __future__ import annotations

import numpy as np
import pytest
from mimosa import EncodedSequences

from motifhorde.evaluation import (
    EvaluatedMotif,
    PerformanceEvaluator,
    _precision_recall_curve,
    _roc_curve,
    _score_cutoff,
    _bootstrap_indices,
    _build_bootstrap_tasks,
)


def test_odd_even_and_seed_order(sequences, tmp_path):
    tasks, test_batches = _build_bootstrap_tasks(
        sequences,
        sequences,
        {"length": [8, 10]},
        str(tmp_path),
        2,
        41,
        [f"p{i}" for i in range(len(sequences))],
        [f"b{i}" for i in range(len(sequences))],
    )

    assert _bootstrap_indices(5, "odd") == ([0, 2, 4], [1, 3])
    assert [task["task_id"] for task in tasks] == [
        "task_0000",
        "task_0001",
        "task_0002",
        "task_0003",
    ]
    assert [task["seed"] for task in tasks] == [42, 43, 44, 45]
    assert [task["split"] for task in tasks] == ["odd", "even", "odd", "even"]
    assert list(test_batches[0][0]) == list(sequences[1])


def test_evaluator_does_not_mutate_mimosa_model(pwm_model, sequences):
    stats = PerformanceEvaluator("peaks").evaluate(
        pwm_model, sequences, sequences, err_threshold=0.1
    )

    assert set(stats) == {"auPRC", "auROC", "pauPRC", "pauROC"}
    assert not hasattr(pwm_model, "config")
    record = EvaluatedMotif(
        pwm_model,
        "task_0000:motif-0001",
        "task_0000",
        "odd",
        {"length": 4},
        stats,
    )
    assert record.params["length"] == 4
    assert record.metrics["auROC"] == stats["auROC"]


def test_curve_reference_values_include_ties_and_inclusive_cutoff():
    labels = np.array([1, 0, 1, 0], dtype=np.int8)
    scores = np.array([0.9, 0.8, 0.8, 0.1], dtype=np.float32)

    tpr, fpr, thresholds = _roc_curve(labels, scores)
    precision, recall, pr_thresholds = _precision_recall_curve(labels, scores)

    np.testing.assert_allclose(tpr, [0.0, 0.5, 1.0, 1.0])
    np.testing.assert_allclose(fpr, [0.0, 0.0, 0.5, 1.0])
    np.testing.assert_allclose(thresholds, [np.inf, 0.9, 0.8, 0.1])
    np.testing.assert_allclose(precision, [1.0, 1.0, 2 / 3, 0.5])
    np.testing.assert_allclose(recall, [0.0, 0.5, 1.0, 1.0])
    np.testing.assert_allclose(pr_thresholds, [np.inf, 0.9, 0.8, 0.1])
    assert _score_cutoff(scores[labels == 0], 0.5) == float(np.float32(0.8))
    assert np.trapezoid(precision, recall) == pytest.approx(11 / 12)
    assert np.trapezoid(tpr, fpr) == pytest.approx(7 / 8)


def test_evaluator_rejects_missing_valid_class_scores(pwm_model, sequences):
    short = EncodedSequences.from_strings(["AC"])
    evaluator = PerformanceEvaluator("peaks")

    with pytest.raises(ValueError, match="positive score"):
        evaluator.evaluate(pwm_model, short, sequences, err_threshold=0.1)
    with pytest.raises(ValueError, match="negative score"):
        evaluator.evaluate(pwm_model, sequences, short, err_threshold=0.1)
