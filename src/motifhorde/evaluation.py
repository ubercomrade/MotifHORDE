"""Performance evaluation and bootstrap validation."""

from __future__ import annotations

import itertools
import logging
import multiprocessing
import os
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from typing import Any, Dict, Iterable, List, Tuple, TypedDict

import numpy as np
from mimosa import GenericModel, scan_model
from mimosa.batches import SequenceBatch, make_sequence_batch, row_values
from mimosa.functions import (
    cut_prc,
    cut_roc,
    format_params,
    lookup_score_for_tail_probability,
    precision_recall_curve,
    roc_curve,
    standardized_pauc,
)
from mimosa.scanning import calculate_threshold_table

from .discovery import MotifDiscoveryTool
from .io import write_fasta
from .log_format import (
    format_elapsed as _format_elapsed,
    format_log_params as _format_log_params,
)

logger = logging.getLogger("evaluation")
bootstrap_logger = logging.getLogger("bootstrap")
bootstrap_discovery_logger = logging.getLogger("bootstrap.discovery")
bootstrap_evaluate_logger = logging.getLogger("bootstrap.evaluate")


class BootstrapTask(TypedDict):
    index: int
    params: dict[str, Any]
    params_suffix: str
    step_name: str
    fg_path: str
    bg_path: str
    output_dir: str
    number_of_motifs: int
    seed: int | None


class BootstrapDiscoveryResult(TypedDict):
    index: int
    params: dict[str, Any]
    params_suffix: str
    step_name: str
    elapsed: float
    motifs: list[GenericModel]


def select_sequence_rows(batch, indices: Iterable[int]):
    """Return a SequenceBatch with selected rows."""
    return make_sequence_batch(row_values(batch, int(index)) for index in indices)


def all_valid_scores(model: GenericModel, sequences) -> np.ndarray:
    """Return all valid best-strand positional scores."""
    scores = scan_model(model, sequences, strand="best")
    return scores["values"][scores["mask"]]


def best_scores(model: GenericModel, sequences) -> np.ndarray:
    """Return the best positional score for each sequence."""
    scores = scan_model(model, sequences, strand="best")
    result = np.full(len(scores["lengths"]), -np.inf, dtype=np.float32)
    for index, length in enumerate(scores["lengths"]):
        if length > 0:
            result[index] = np.max(scores["values"][index, : int(length)])
    return result[np.isfinite(result)]


class PerformanceEvaluator:
    """Compute binary classification metrics for motif models."""

    def __init__(self, background_type: str = "sites") -> None:
        self.background_type = background_type

    def evaluate(
        self,
        motif: GenericModel,
        positives,
        negatives,
        err_threshold: float,
    ) -> Dict[str, Any]:
        true_max_scores = best_scores(motif, positives)

        if self.background_type == "sites":
            false_scores = all_valid_scores(motif, negatives)
        elif self.background_type == "peaks":
            false_scores = best_scores(motif, negatives)
        else:
            logger.warning(
                "invalid_background_type | value=%s | fallback=peaks",
                self.background_type,
            )
            false_scores = best_scores(motif, negatives)

        classification = np.concatenate(
            (
                np.ones(len(true_max_scores), dtype=np.int8),
                np.zeros(len(false_scores), dtype=np.int8),
            )
        )
        scores = np.concatenate(
            (
                true_max_scores.astype(np.float32, copy=False),
                false_scores.astype(np.float32, copy=False),
            )
        )

        prec, rec, uniq_scores_pr = precision_recall_curve(classification, scores)
        tpr, fpr, uniq_scores_roc = roc_curve(classification, scores)

        auprc = float(np.trapezoid(prec, rec))
        auroc = float(np.trapezoid(tpr, fpr))

        threshold_table = calculate_threshold_table(motif, negatives, strand="best")
        score_cutoff = lookup_score_for_tail_probability(threshold_table, err_threshold)

        tpr_cut, fpr_cut, _ = cut_roc(tpr, fpr, uniq_scores_roc, score_cutoff)
        pauroc_raw = float(np.trapezoid(tpr_cut, fpr_cut))

        rec_cut, prec_cut, _ = cut_prc(rec, prec, uniq_scores_pr, score_cutoff)
        pauprc_raw = float(np.trapezoid(prec_cut, rec_cut))

        e = float(fpr_cut[-1]) if len(fpr_cut) else 0.0
        r = float(rec_cut[-1]) if len(rec_cut) else 0.0
        pauroc = standardized_pauc(pauroc_raw, pauc_min=(e * e / 2.0), pauc_max=e)
        pauprc = standardized_pauc(pauprc_raw, pauc_min=(0.5 * r), pauc_max=r)

        stats = {
            "auPRC": auprc,
            "auROC": auroc,
            "pauPRC": pauprc,
            "pauROC": pauroc,
        }
        motif.config["statistics"] = stats

        return {
            "PRC": {"RECALL": rec.tolist(), "PRECISION": prec.tolist()},
            "ROC": {"FPR": fpr.tolist(), "TPR": tpr.tolist()},
            **stats,
        }


class Bootstrapper:
    """Run odd/even bootstrap discovery and evaluation."""

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
        peaks,
        background,
        number_of_motifs: int,
        err_threshold: float,
        discovery_params: Dict[str, Iterable[Any]],
    ) -> Tuple[Dict[str, Any], List[GenericModel]]:
        start_time = time.perf_counter()
        task_parent = os.path.join(self.output_dir, self.discovery_tool.name)
        os.makedirs(task_parent, exist_ok=True)
        with tempfile.TemporaryDirectory(
            dir=task_parent,
            prefix="bootstrap_parallel_",
        ) as task_root:
            tasks, test_batches = _build_bootstrap_tasks(
                peaks,
                background,
                discovery_params,
                task_root,
                number_of_motifs,
                self.seed,
            )
            bootstrap_logger.info(
                "start | param_sets=%d | tasks=%d | jobs=%d | motifs_per_task=%d | fpr=%s",
                len(tasks) // 2,
                len(tasks),
                self.jobs,
                number_of_motifs,
                err_threshold,
            )
            discovery_start = time.perf_counter()
            results = _run_bootstrap_discovery_tasks(
                self.discovery_tool,
                tasks,
                self.jobs,
            )
            sorted_results = sorted(results, key=lambda item: item["index"])
            motifs_found = sum(len(result["motifs"]) for result in sorted_results)
            for task_number, result in enumerate(sorted_results, start=1):
                bootstrap_discovery_logger.info(
                    "done | task=%d/%d | split=%s | params=%s | motifs_found=%d | elapsed=%s",
                    task_number,
                    len(sorted_results),
                    result["step_name"],
                    _format_log_params(result["params"]),
                    len(result["motifs"]),
                    _format_elapsed(result["elapsed"]),
                )
            bootstrap_discovery_logger.info(
                "complete | tasks=%d | motifs_found=%d | elapsed=%s",
                len(sorted_results),
                motifs_found,
                _format_elapsed(time.perf_counter() - discovery_start),
            )

            statistics, bootstrap_motifs = _evaluate_bootstrap_results(
                self.evaluator,
                sorted_results,
                test_batches,
                background,
                err_threshold,
            )
            bootstrap_evaluate_logger.info(
                "done | motifs=%d | statistics=%d",
                len(bootstrap_motifs),
                len(statistics),
            )
            bootstrap_logger.info(
                "done | elapsed=%s",
                _format_elapsed(time.perf_counter() - start_time),
            )
            return statistics, bootstrap_motifs


def _bootstrap_indices(n_peaks: int, step_name: str) -> tuple[list[int], list[int]]:
    if step_name == "odd":
        train_indices = [index for index in range(n_peaks) if (index + 1) % 2 != 0]
        test_indices = [index for index in range(n_peaks) if (index + 1) % 2 == 0]
        return train_indices, test_indices
    if step_name == "even":
        train_indices = [index for index in range(n_peaks) if (index + 1) % 2 == 0]
        test_indices = [index for index in range(n_peaks) if (index + 1) % 2 != 0]
        return train_indices, test_indices
    raise ValueError(f"Unknown bootstrap split: {step_name}")


def _task_seed(base_seed: int | None, task_index: int) -> int | None:
    if base_seed is None:
        return None
    return base_seed + task_index + 1


def _safe_task_component(value: str) -> str:
    return value.replace("/", "_").replace("\\", "_").replace(":", "-")


def _build_bootstrap_tasks(
    peaks: SequenceBatch,
    background: SequenceBatch,
    discovery_params: Dict[str, Iterable[Any]],
    task_root: str,
    number_of_motifs: int,
    seed: int | None,
) -> tuple[list[BootstrapTask], dict[int, SequenceBatch]]:
    tasks: list[BootstrapTask] = []
    test_batches: dict[int, SequenceBatch] = {}
    n_peaks = len(peaks["lengths"])
    param_keys = list(discovery_params.keys())
    param_values = list(discovery_params.values())

    for combination in itertools.product(*param_values):
        current_params = dict(zip(param_keys, combination))
        params_suffix = format_params(current_params)

        for step_name in ["odd", "even"]:
            task_index = len(tasks)
            train_indices, test_indices = _bootstrap_indices(n_peaks, step_name)
            train_peaks = select_sequence_rows(peaks, train_indices)
            test_batches[task_index] = select_sequence_rows(peaks, test_indices)

            task_dir = os.path.join(
                task_root,
                f"task_{task_index:04d}_{_safe_task_component(params_suffix)}_{step_name}",
            )
            os.makedirs(task_dir, exist_ok=True)
            fg_path = os.path.join(task_dir, "train.fasta")
            bg_path = os.path.join(task_dir, "background.fasta")
            write_fasta(train_peaks, fg_path)
            write_fasta(background, bg_path)

            tasks.append(
                {
                    "index": task_index,
                    "params": current_params,
                    "params_suffix": params_suffix,
                    "step_name": step_name,
                    "fg_path": fg_path,
                    "bg_path": bg_path,
                    "output_dir": task_dir,
                    "number_of_motifs": number_of_motifs,
                    "seed": _task_seed(seed, task_index),
                }
            )

    return tasks, test_batches


def _run_discovery_task(
    discovery_tool: MotifDiscoveryTool,
    task: BootstrapTask,
) -> BootstrapDiscoveryResult:
    start_time = time.perf_counter()
    try:
        kwargs = dict(task["params"])
        if task["seed"] is not None and discovery_tool.name == "sitega":
            kwargs["seed"] = task["seed"]

        motifs = discovery_tool.discover(
            task["fg_path"],
            task["bg_path"],
            task["output_dir"],
            number_of_motifs=task["number_of_motifs"],
            **kwargs,
        )
    except Exception as exc:
        raise RuntimeError(
            f"Bootstrap discovery failed for params={task['params']} "
            f"split={task['step_name']}"
        ) from exc

    return {
        "index": task["index"],
        "params": task["params"],
        "params_suffix": task["params_suffix"],
        "step_name": task["step_name"],
        "elapsed": time.perf_counter() - start_time,
        "motifs": motifs,
    }


def _run_bootstrap_discovery_tasks(
    discovery_tool: MotifDiscoveryTool,
    tasks: list[BootstrapTask],
    jobs: int,
) -> list[BootstrapDiscoveryResult]:
    if jobs == 1:
        return [_run_discovery_task(discovery_tool, task) for task in tasks]

    mp_context = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(max_workers=jobs, mp_context=mp_context) as executor:
        return list(executor.map(_run_discovery_task, repeat(discovery_tool), tasks))


def _evaluate_bootstrap_results(
    evaluator: PerformanceEvaluator,
    results: list[BootstrapDiscoveryResult],
    test_batches: dict[int, SequenceBatch],
    background: SequenceBatch,
    err_threshold: float,
) -> tuple[Dict[str, Any], List[GenericModel]]:
    statistics: Dict[str, Any] = {}
    bootstrap_motifs: List[GenericModel] = []

    for result in sorted(results, key=lambda item: item["index"]):
        test_peaks = test_batches[result["index"]]
        for motif in result["motifs"]:
            stats = evaluator.evaluate(motif, test_peaks, background, err_threshold)
            motif.name = f"{motif.name}_{result['params_suffix']}_{result['step_name']}"
            statistics[motif.name] = stats
            bootstrap_motifs.append(motif)

    return statistics, bootstrap_motifs
