"""Motif evaluation and deterministic odd/even bootstrap discovery."""

from __future__ import annotations

import itertools
import multiprocessing
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from itertools import repeat
from types import MappingProxyType
from typing import Any, Iterable, TypedDict

import numpy as np
from mimosa import EncodedSequences, MotifModel, scan

from .discovery import MotifDiscoveryTool
from .io import write_fasta

VALIDATION_METRICS = ("auPRC", "auROC", "pauPRC", "pauROC")


@dataclass(frozen=True, slots=True)
class EvaluatedMotif:
    """A discovered model plus bootstrap context and validation metrics."""

    model: MotifModel
    id: str
    task_id: str
    split: str
    params: MappingProxyType
    metrics: MappingProxyType

    def __post_init__(self) -> None:
        object.__setattr__(self, "params", MappingProxyType(dict(self.params)))
        object.__setattr__(
            self,
            "metrics",
            MappingProxyType(
                {key: float(value) for key, value in self.metrics.items()}
            ),
        )


class BootstrapTask(TypedDict):
    task_id: str
    index: int
    params: dict[str, Any]
    split: str
    fg_path: str
    bg_path: str
    output_dir: str
    number_of_motifs: int
    seed: int | None


def select_sequence_rows(
    batch: EncodedSequences, indices: Iterable[int]
) -> EncodedSequences:
    """Select rows without changing their order."""
    return EncodedSequences.from_rows(batch[int(index)] for index in indices)


def _valid_positions(model: MotifModel, sequences: EncodedSequences) -> np.ndarray:
    return scan(model, sequences, strands="best").data


def all_valid_scores(model: MotifModel, sequences: EncodedSequences) -> np.ndarray:
    """Return all finite positional best-strand scores."""
    scores = _valid_positions(model, sequences)
    return scores[np.isfinite(scores)]


def best_scores(model: MotifModel, sequences: EncodedSequences) -> np.ndarray:
    """Return the best finite positional score for every sequence."""
    positions = scan(model, sequences, strands="best")
    values: list[float] = []
    for index in range(len(sequences)):
        row = positions[index]
        if row.size:
            best = float(np.max(row))
            if np.isfinite(best):
                values.append(best)
    return np.asarray(values, dtype=np.float32)


def _roc_curve(labels: np.ndarray, scores: np.ndarray):
    if len(scores) == 0:
        return np.array([0.0]), np.array([0.0]), np.array([np.inf])
    indexes = np.argsort(-scores, kind="stable")
    sorted_scores = scores[indexes]
    sorted_labels = labels[indexes]
    size = np.unique(scores).size + 1
    tpr = np.zeros(size)
    fpr = np.zeros(size)
    thresholds = np.zeros(size)
    thresholds[0] = np.inf
    true_count = false_count = 0
    positives = np.sum(labels == 1)
    negatives = np.sum(labels == 0)
    position = 1
    current_score = sorted_scores[0]
    for index, score in enumerate(sorted_scores):
        if sorted_labels[index] == 1:
            true_count += 1
        else:
            false_count += 1
        if index == len(scores) - 1 or current_score != sorted_scores[index + 1]:
            thresholds[position] = score
            tpr[position] = true_count / positives if positives else 0.0
            fpr[position] = false_count / negatives if negatives else 0.0
            position += 1
            if index < len(scores) - 1:
                current_score = sorted_scores[index + 1]
    return tpr[:position], fpr[:position], thresholds[:position]


def _precision_recall_curve(labels: np.ndarray, scores: np.ndarray):
    if len(scores) == 0:
        return np.array([1.0]), np.array([0.0]), np.array([np.inf])
    indexes = np.argsort(-scores, kind="stable")
    sorted_scores = scores[indexes]
    sorted_labels = labels[indexes]
    size = np.unique(scores).size + 1
    precision = np.zeros(size)
    recall = np.zeros(size)
    thresholds = np.zeros(size)
    precision[0] = 1.0
    thresholds[0] = np.inf
    true_count = false_count = 0
    positives = np.sum(labels == 1)
    negatives = np.sum(labels == 0)
    ratio = 1.0 if negatives == 0 else positives / negatives
    position = 1
    current_score = sorted_scores[0]
    for index, score in enumerate(sorted_scores):
        if sorted_labels[index] == 1:
            true_count += 1
        else:
            false_count += 1
        if index == len(scores) - 1 or current_score != sorted_scores[index + 1]:
            thresholds[position] = score
            denominator = true_count + ratio * false_count
            precision[position] = true_count / denominator if denominator else 1.0
            recall[position] = true_count / positives if positives else 0.0
            position += 1
            if index < len(scores) - 1:
                current_score = sorted_scores[index + 1]
    return precision[:position], recall[:position], thresholds[:position]


def _cut_roc(first, second, thresholds, cutoff):
    if cutoff == -np.inf:
        return first, second, thresholds
    mask = thresholds >= cutoff
    if not np.any(mask):
        return (
            np.array([first[0]], dtype=first.dtype),
            np.array([0.0], dtype=second.dtype),
            np.array([cutoff], dtype=thresholds.dtype),
        )
    last = int(np.where(mask)[0][-1])
    if thresholds[last] == cutoff or last == len(thresholds) - 1:
        return first[: last + 1], second[: last + 1], thresholds[: last + 1]
    s0, s1 = float(thresholds[last]), float(thresholds[last + 1])
    alpha = 0.0 if s0 == s1 else (cutoff - s0) / (s1 - s0)
    first_cut = float(first[last]) + alpha * (
        float(first[last + 1]) - float(first[last])
    )
    second_cut = float(second[last]) + alpha * (
        float(second[last + 1]) - float(second[last])
    )
    return (
        np.concatenate([first[: last + 1], np.array([first_cut], dtype=first.dtype)]),
        np.concatenate(
            [second[: last + 1], np.array([second_cut], dtype=second.dtype)]
        ),
        np.concatenate(
            [thresholds[: last + 1], np.array([cutoff], dtype=thresholds.dtype)]
        ),
    )


def _cut_prc(recall, precision, thresholds, cutoff):
    if cutoff == -np.inf:
        return recall, precision, thresholds
    mask = thresholds >= cutoff
    if not np.any(mask):
        return (
            np.array([0.0], dtype=recall.dtype),
            np.array([precision[0]], dtype=precision.dtype),
            np.array([cutoff], dtype=thresholds.dtype),
        )
    last = int(np.where(mask)[0][-1])
    if thresholds[last] == cutoff or last == len(thresholds) - 1:
        return recall[: last + 1], precision[: last + 1], thresholds[: last + 1]
    s0, s1 = float(thresholds[last]), float(thresholds[last + 1])
    alpha = 0.0 if s0 == s1 else (cutoff - s0) / (s1 - s0)
    recall_cut = float(recall[last]) + alpha * (
        float(recall[last + 1]) - float(recall[last])
    )
    precision_cut = float(precision[last]) + alpha * (
        float(precision[last + 1]) - float(precision[last])
    )
    return (
        np.concatenate(
            [recall[: last + 1], np.array([recall_cut], dtype=recall.dtype)]
        ),
        np.concatenate(
            [precision[: last + 1], np.array([precision_cut], dtype=precision.dtype)]
        ),
        np.concatenate(
            [thresholds[: last + 1], np.array([cutoff], dtype=thresholds.dtype)]
        ),
    )


def _score_cutoff(scores: np.ndarray, tail_probability: float) -> float:
    if scores.size == 0:
        raise ValueError("at least one negative positional score is required")
    if (
        not np.isfinite(tail_probability)
        or tail_probability < 0
        or tail_probability > 1
    ):
        raise ValueError("false-positive threshold must lie in [0, 1]")
    sorted_scores = np.sort(np.asarray(scores, dtype=np.float32))[::-1]
    unique_scores, counts = np.unique(sorted_scores, return_counts=True)
    unique_scores = unique_scores[::-1]
    counts = counts[::-1]
    log_tail = -np.log10(np.cumsum(counts) / scores.size)
    if tail_probability <= 0:
        return float(unique_scores[0])
    mask = log_tail >= -np.log10(tail_probability)
    return (
        float(unique_scores[np.where(mask)[0][-1]])
        if np.any(mask)
        else float(unique_scores[-1])
    )


def _standardized_pauc(raw: float, minimum: float, maximum: float) -> float:
    denominator = maximum - minimum
    return 0.5 if denominator <= 0 else 0.5 * (1.0 + (raw - minimum) / denominator)


class PerformanceEvaluator:
    """Compute ROC, PR, and standardized partial-area metrics."""

    def __init__(self, background_type: str = "sites") -> None:
        if background_type not in {"sites", "peaks"}:
            raise ValueError("background_type must be 'sites' or 'peaks'")
        self.background_type = background_type

    def evaluate(
        self,
        motif: MotifModel,
        positives: EncodedSequences,
        negatives: EncodedSequences,
        err_threshold: float,
    ) -> dict[str, float]:
        if not np.isfinite(err_threshold):
            raise ValueError("false-positive threshold must be finite")
        positive_scores = best_scores(motif, positives)
        negative_scores = (
            all_valid_scores(motif, negatives)
            if self.background_type == "sites"
            else best_scores(motif, negatives)
        )
        if positive_scores.size == 0:
            raise ValueError("motif evaluation requires at least one positive score")
        if negative_scores.size == 0:
            raise ValueError("motif evaluation requires at least one negative score")

        labels = np.concatenate(
            [
                np.ones(positive_scores.size, dtype=np.int8),
                np.zeros(negative_scores.size, dtype=np.int8),
            ]
        )
        scores = np.concatenate([positive_scores, negative_scores]).astype(
            np.float32, copy=False
        )
        precision, recall, pr_thresholds = _precision_recall_curve(labels, scores)
        tpr, fpr, roc_thresholds = _roc_curve(labels, scores)
        cutoff = _score_cutoff(negative_scores, float(err_threshold))
        tpr_cut, fpr_cut, _ = _cut_roc(tpr, fpr, roc_thresholds, cutoff)
        recall_cut, precision_cut, _ = _cut_prc(
            recall, precision, pr_thresholds, cutoff
        )
        fpr_end = float(fpr_cut[-1]) if fpr_cut.size else 0.0
        recall_end = float(recall_cut[-1]) if recall_cut.size else 0.0
        return {
            "auPRC": float(np.trapezoid(precision, recall)),
            "auROC": float(np.trapezoid(tpr, fpr)),
            "pauPRC": _standardized_pauc(
                float(np.trapezoid(precision_cut, recall_cut)),
                0.5 * recall_end,
                recall_end,
            ),
            "pauROC": _standardized_pauc(
                float(np.trapezoid(tpr_cut, fpr_cut)),
                fpr_end * fpr_end / 2.0,
                fpr_end,
            ),
        }


class Bootstrapper:
    """Run deterministic odd/even discovery tasks."""

    def __init__(
        self,
        discovery_tool: MotifDiscoveryTool,
        evaluator: PerformanceEvaluator,
        output_dir: str,
        jobs: int = 1,
        seed: int | None = None,
    ) -> None:
        self.discovery_tool = discovery_tool
        self.evaluator = evaluator
        self.output_dir = output_dir
        self.jobs = jobs
        self.seed = seed

    def run(
        self,
        peaks: EncodedSequences,
        background: EncodedSequences,
        number_of_motifs: int,
        err_threshold: float,
        discovery_params: dict[str, Iterable[Any]],
        peak_names: Iterable[str] | None = None,
        background_names: Iterable[str] | None = None,
    ) -> list[EvaluatedMotif]:
        peak_names = tuple(peak_names or (f"peak-{i}" for i in range(len(peaks))))
        background_names = tuple(
            background_names or (f"background-{i}" for i in range(len(background)))
        )
        task_root = os.path.join(
            self.output_dir, self.discovery_tool.name, "bootstrap", "tasks"
        )
        os.makedirs(task_root, exist_ok=True)
        tasks, test_batches = _build_bootstrap_tasks(
            peaks,
            background,
            discovery_params,
            task_root,
            number_of_motifs,
            self.seed,
            peak_names,
            background_names,
        )
        results = _run_bootstrap_discovery_tasks(self.discovery_tool, tasks, self.jobs)
        return _evaluate_bootstrap_results(
            self.evaluator, results, test_batches, background, err_threshold
        )


def _bootstrap_indices(n_peaks: int, split: str) -> tuple[list[int], list[int]]:
    if split == "odd":
        train = [i for i in range(n_peaks) if (i + 1) % 2]
        return train, [i for i in range(n_peaks) if not (i + 1) % 2]
    if split == "even":
        train = [i for i in range(n_peaks) if not (i + 1) % 2]
        return train, [i for i in range(n_peaks) if (i + 1) % 2]
    raise ValueError(f"Unknown bootstrap split: {split}")


def _task_seed(base_seed: int | None, task_index: int) -> int | None:
    return None if base_seed is None else base_seed + task_index + 1


def _build_bootstrap_tasks(
    peaks: EncodedSequences,
    background: EncodedSequences,
    discovery_params: dict[str, Iterable[Any]],
    task_root: str,
    number_of_motifs: int,
    seed: int | None,
    peak_names: Iterable[str] | None = None,
    background_names: Iterable[str] | None = None,
) -> tuple[list[BootstrapTask], dict[int, EncodedSequences]]:
    peak_names = tuple(peak_names or (f"peak-{i}" for i in range(len(peaks))))
    background_names = tuple(
        background_names or (f"background-{i}" for i in range(len(background)))
    )
    if len(peak_names) != len(peaks) or len(background_names) != len(background):
        raise ValueError("FASTA name counts must match sequence counts")
    keys = list(discovery_params)
    tasks: list[BootstrapTask] = []
    test_batches: dict[int, EncodedSequences] = {}
    for values in itertools.product(*(discovery_params[key] for key in keys)):
        params = dict(zip(keys, values, strict=True))
        for split in ("odd", "even"):
            index = len(tasks)
            train_indices, test_indices = _bootstrap_indices(len(peaks), split)
            task_id = f"task_{index:04d}"
            task_dir = os.path.join(task_root, task_id)
            os.makedirs(task_dir, exist_ok=True)
            train_peaks = select_sequence_rows(peaks, train_indices)
            test_batches[index] = select_sequence_rows(peaks, test_indices)
            fg_path = os.path.join(task_dir, "train.fasta")
            bg_path = os.path.join(task_dir, "background.fasta")
            write_fasta(train_peaks, fg_path, [peak_names[i] for i in train_indices])
            write_fasta(background, bg_path, background_names)
            tasks.append(
                {
                    "task_id": task_id,
                    "index": index,
                    "params": params,
                    "split": split,
                    "fg_path": fg_path,
                    "bg_path": bg_path,
                    "output_dir": task_dir,
                    "number_of_motifs": number_of_motifs,
                    "seed": _task_seed(seed, index),
                }
            )
    return tasks, test_batches


def _run_discovery_task(
    discovery_tool: MotifDiscoveryTool, task: BootstrapTask
) -> dict[str, Any]:
    kwargs = dict(task["params"])
    if task["seed"] is not None and discovery_tool.name == "sitega":
        kwargs["seed"] = task["seed"]
    try:
        motifs = discovery_tool.discover(
            task["fg_path"],
            task["bg_path"],
            task["output_dir"],
            number_of_motifs=task["number_of_motifs"],
            **kwargs,
        )
    except Exception as exc:
        raise RuntimeError(
            f"Bootstrap discovery failed for task={task['task_id']} "
            f"split={task['split']} params={task['params']}"
        ) from exc
    if not motifs:
        raise RuntimeError(
            f"Bootstrap discovery returned no models for {task['task_id']}"
        )
    return {"task": task, "motifs": motifs}


def _run_bootstrap_discovery_tasks(
    discovery_tool: MotifDiscoveryTool, tasks: list[BootstrapTask], jobs: int
) -> list[dict[str, Any]]:
    if jobs < 1:
        raise ValueError("jobs must be positive")
    if jobs == 1:
        return [_run_discovery_task(discovery_tool, task) for task in tasks]
    context = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(max_workers=jobs, mp_context=context) as executor:
        return list(executor.map(_run_discovery_task, repeat(discovery_tool), tasks))


def _evaluate_bootstrap_results(
    evaluator: PerformanceEvaluator,
    results: list[dict[str, Any]],
    test_batches: dict[int, EncodedSequences],
    background: EncodedSequences,
    err_threshold: float,
) -> list[EvaluatedMotif]:
    evaluated: list[EvaluatedMotif] = []
    for result in sorted(results, key=lambda item: item["task"]["index"]):
        task = result["task"]
        for index, model in enumerate(result["motifs"], start=1):
            metrics = evaluator.evaluate(
                model, test_batches[task["index"]], background, err_threshold
            )
            evaluated.append(
                EvaluatedMotif(
                    model=model,
                    id=f"{task['task_id']}:motif-{index:04d}",
                    task_id=task["task_id"],
                    split=task["split"],
                    params=task["params"],
                    metrics=metrics,
                )
            )
    return evaluated
