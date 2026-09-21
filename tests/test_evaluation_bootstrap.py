from __future__ import annotations

from motifhorde.evaluation import (
    EvaluatedMotif,
    PerformanceEvaluator,
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
