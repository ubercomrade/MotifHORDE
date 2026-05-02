"""Performance evaluation and bootstrap validation."""

from __future__ import annotations

import itertools
import os
import tempfile
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np

from .batches import make_sequence_batch, row_values
from .discovery import MotifDiscoveryTool
from .functions import (
    cut_prc,
    cut_roc,
    format_params,
    lookup_score_for_tail_probability,
    precision_recall_curve,
    roc_curve,
    standardized_pauc,
)
from .io import write_fasta
from .models import GenericModel, calculate_threshold_table, scan_model


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
            print(f"Incorrect background_type: {self.background_type}, set as `peaks`")
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

        auprc = float(np.trapz(prec, rec))
        auroc = float(np.trapz(tpr, fpr))

        threshold_table = calculate_threshold_table(motif, negatives, strand="best")
        score_cutoff = lookup_score_for_tail_probability(threshold_table, err_threshold)

        tpr_cut, fpr_cut, _ = cut_roc(tpr, fpr, uniq_scores_roc, score_cutoff)
        pauroc_raw = float(np.trapz(tpr_cut, fpr_cut))

        rec_cut, prec_cut, _ = cut_prc(rec, prec, uniq_scores_pr, score_cutoff)
        pauprc_raw = float(np.trapz(prec_cut, rec_cut))

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

    def __init__(self, discovery_tool: MotifDiscoveryTool, evaluator: PerformanceEvaluator, output_dir: str) -> None:
        self.discovery_tool = discovery_tool
        self.evaluator = evaluator
        self.output_dir = output_dir

    def run(
        self,
        peaks,
        background,
        number_of_motifs: int,
        err_threshold: float,
        discovery_params: Dict[str, Iterable[Any]],
    ) -> Tuple[Dict[str, Any], List[GenericModel]]:
        statistics: Dict[str, Any] = {}
        bootstrap_motifs: List[GenericModel] = []

        print("Starting bootstrap...")

        param_keys = list(discovery_params.keys())
        param_values = list(discovery_params.values())

        for combination in itertools.product(*param_values):
            current_params = dict(zip(param_keys, combination))
            params_suffix = format_params(current_params)

            for step_name in ["odd", "even"]:
                prefix = f"bootstrap_{params_suffix}_{step_name}"
                with tempfile.TemporaryDirectory(
                    dir=os.path.join(self.output_dir, self.discovery_tool.name),
                    prefix=prefix,
                    delete=True,
                ) as tmp_dir:
                    n_peaks = len(peaks["lengths"])
                    if step_name == "odd":
                        train_indices = [i for i in range(n_peaks) if (i + 1) % 2 != 0]
                        test_indices = [i for i in range(n_peaks) if (i + 1) % 2 == 0]
                    else:
                        train_indices = [i for i in range(n_peaks) if (i + 1) % 2 == 0]
                        test_indices = [i for i in range(n_peaks) if (i + 1) % 2 != 0]

                    train_peaks = select_sequence_rows(peaks, train_indices)
                    test_peaks = select_sequence_rows(peaks, test_indices)

                    fg_path = os.path.join(tmp_dir, "train.fasta")
                    bg_path = os.path.join(tmp_dir, "background.fasta")
                    write_fasta(train_peaks, fg_path)
                    write_fasta(background, bg_path)

                    motifs = self.discovery_tool.discover(
                        fg_path,
                        bg_path,
                        tmp_dir,
                        number_of_motifs=number_of_motifs,
                        **current_params,
                    )

                    for motif in motifs:
                        stats = self.evaluator.evaluate(motif, test_peaks, background, err_threshold)
                        motif.name = f"{motif.name}_{params_suffix}_{step_name}"
                        statistics[motif.name] = stats

                    bootstrap_motifs.extend(motifs)

        return statistics, bootstrap_motifs
