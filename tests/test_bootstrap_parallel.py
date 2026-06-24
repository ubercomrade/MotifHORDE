from __future__ import annotations

import os

import numpy as np
from mimosa.batches import make_sequence_batch

from motifhorde.evaluation import (
    Bootstrapper,
    _bootstrap_indices,
    _build_bootstrap_tasks,
    _task_seed,
)
from motifhorde.models import GenericModel


def _batch(*rows: list[int]):
    return make_sequence_batch([np.asarray(row, dtype=np.int8) for row in rows])


class PickleableDiscoveryTool:
    name = "fake"

    def discover(
        self,
        foreground: str,
        background: str,
        output_dir: str,
        number_of_motifs: int,
        **kwargs,
    ):
        assert os.path.exists(foreground)
        assert os.path.exists(background)
        return [
            GenericModel(
                "test",
                f"Fake-{kwargs['length']}-{index + 1}",
                None,
                int(kwargs["length"]),
                {},
            )
            for index in range(number_of_motifs)
        ]


class DeterministicEvaluator:
    def evaluate(self, motif, positives, negatives, err_threshold):
        value = float(motif.length + len(positives["lengths"]))
        stats = {
            "auPRC": value,
            "auROC": value,
            "pauPRC": value,
            "pauROC": value,
        }
        motif.config["statistics"] = stats
        return stats


def test_bootstrap_indices_match_odd_even_splits():
    assert _bootstrap_indices(5, "odd") == ([0, 2, 4], [1, 3])
    assert _bootstrap_indices(5, "even") == ([1, 3], [0, 2, 4])


def test_build_bootstrap_tasks_has_stable_order_and_unique_dirs(tmp_path):
    peaks = _batch([0, 1, 2], [1, 2, 3], [2, 3, 0], [3, 0, 1])
    background = _batch([0, 0, 0], [1, 1, 1])

    tasks, test_batches = _build_bootstrap_tasks(
        peaks,
        background,
        {"length": [8, 10]},
        os.fspath(tmp_path),
        number_of_motifs=2,
        seed=100,
    )

    assert [task["index"] for task in tasks] == [0, 1, 2, 3]
    assert [task["step_name"] for task in tasks] == ["odd", "even", "odd", "even"]
    assert [task["params"] for task in tasks] == [
        {"length": 8},
        {"length": 8},
        {"length": 10},
        {"length": 10},
    ]
    assert [task["seed"] for task in tasks] == [101, 102, 103, 104]
    assert len({task["output_dir"] for task in tasks}) == 4
    assert all(os.path.exists(task["fg_path"]) for task in tasks)
    assert all(os.path.exists(task["bg_path"]) for task in tasks)
    assert set(test_batches) == {0, 1, 2, 3}


def test_task_seed_preserves_nondeterministic_mode_without_base_seed():
    assert _task_seed(None, 0) is None
    assert _task_seed(42, 0) == 43
    assert _task_seed(42, 2) == 45


def test_bootstrap_sequential_and_process_modes_match(tmp_path):
    peaks = _batch([0, 1, 2], [1, 2, 3], [2, 3, 0], [3, 0, 1])
    background = _batch([0, 0, 0], [1, 1, 1])
    discovery_params = {"length": [8, 10]}

    sequential = Bootstrapper(
        PickleableDiscoveryTool(),
        DeterministicEvaluator(),
        os.fspath(tmp_path / "seq"),
        jobs=1,
    ).run(
        peaks,
        background,
        number_of_motifs=1,
        err_threshold=0.001,
        discovery_params=discovery_params,
    )
    parallel = Bootstrapper(
        PickleableDiscoveryTool(),
        DeterministicEvaluator(),
        os.fspath(tmp_path / "par"),
        jobs=2,
    ).run(
        peaks,
        background,
        number_of_motifs=1,
        err_threshold=0.001,
        discovery_params=discovery_params,
    )

    sequential_stats, sequential_motifs = sequential
    parallel_stats, parallel_motifs = parallel
    assert [motif.name for motif in sequential_motifs] == [
        "Fake-8-1_length-8_odd",
        "Fake-8-1_length-8_even",
        "Fake-10-1_length-10_odd",
        "Fake-10-1_length-10_even",
    ]
    assert [motif.name for motif in parallel_motifs] == [
        motif.name for motif in sequential_motifs
    ]
    assert parallel_stats == sequential_stats
