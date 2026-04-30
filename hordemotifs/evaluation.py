"""
evaluation
==========

Classes for computing performance metrics of motif models.  The
evaluation measures implemented here focus on binary classification
between foreground and background sequences using area under the ROC
curve (AUROC), area under the precision–recall curve (AUPRC) and
their partial variants.  These metrics are widely used to assess
motif quality and have different sensitivity to class imbalance.
"""

from __future__ import annotations
from typing import Dict, List, Tuple, Any, Iterable
import os
import tempfile
import itertools
import numpy as np
from .discovery import MotifDiscoveryTool
from .models import MotifModel
from .ragged import RaggedData
from .io import write_fasta
from .functions import (
    precision_recall_curve,
    roc_curve,
    format_params,
    standardized_pauc,
    cut_roc, cut_prc
)
from .ragged import ragged_from_list

class PerformanceEvaluator:
    """Evaluator for motifs.

    Uses the optimised scoring kernel from the original pipeline to
    compute per‑sequence scores and then derives precision–recall and
    ROC curves.  Partial AUC is computed up to a specified FPR
    threshold.
    """

    def __init__(self, background_type: str = 'sites') -> None:
        self.background_type = background_type

    def evaluate(self, motif: MotifModel, positives: RaggedData, negatives: RaggedData,
                 err_threshold: float) -> Dict[str, Any]:
        """Compute performance metrics for a motif model using binary classification.

        This method evaluates motif performance by treating motif scanning scores as predictions
        for binary classification between foreground (positive) and background (negative) sequences.
        Standard metrics including AUROC, AUPRC and their partial variants are computed.

        Parameters
        ----------
        motif : MotifModel
            The motif model to be evaluated.
        positives : RaggedData
            Foreground sequences representing true binding events.
        negatives : RaggedData
            Background sequences representing non-binding events.
        err_threshold : float
            Expectation Recognition Rate (ERR/FPR) threshold used for partial AUC computation.

        Returns
        -------
        Dict[str, Any]
            Dictionary containing performance metrics including:
            - PRC: Precision-Recall curve data
            - ROC: ROC curve data
            - auPRC: Area Under the Precision-Recall Curve
            - auROC: Area Under the ROC Curve
            - pauPRC: Partial Area Under the Precision-Recall Curve
            - pauROC: Partial Area Under the ROC Curve
        """
        # Batch scoring using RaggedData API
        # best_scores returns np.ndarray of shape (n_seq,)
        true_max_scores = motif.best_scores(positives, strand="best")
        
        if self.background_type == 'peaks':
            false_max_scores = motif.best_scores(negatives, strand="best")
        elif self.background_type == 'sites':
            # scan returns RaggedData with scores for ALL positions (reduced to 1D via "best" strand)
            scan_results = motif.scan(negatives, strand="best")
            false_max_scores = scan_results.data # All positions scores
        else:
            print(f'Incorrect background_type: {self.background_type}, set as `peaks`')
            false_max_scores = motif.best_scores(negatives, strand="best")

        # Classification labels: 1 for true positives, 0 for false positives
        classification = np.concatenate((
            np.ones(len(true_max_scores), dtype=np.int8),
            np.zeros(len(false_max_scores), dtype=np.int8)
        ))
        scores = np.concatenate((
            true_max_scores.astype(np.float32),
            false_max_scores.astype(np.float32)
        ))

        # Compute full PRC
        prec, rec, uniq_scores_pr = precision_recall_curve(classification, scores)

        # Compute full ROC
        tpr, fpr, uniq_scores_roc = roc_curve(classification, scores)

        # Full AUC values
        auprc = float(np.trapz(prec, rec))
        auroc = float(np.trapz(tpr, fpr))

        # Compute score cutoff for partial AUC calculations
        score_cutoff = motif._frequency_to_score(err_threshold, background_data=negatives)

        tpr_cut, fpr_cut, thr_roc_cut = cut_roc(tpr, fpr, uniq_scores_roc, score_cutoff)
        pauroc_raw = float(np.trapz(tpr_cut, fpr_cut))

        rec_cut, prec_cut, thr_pr_cut = cut_prc(rec, prec, uniq_scores_pr, score_cutoff)
        pauprc_raw = float(np.trapz(prec_cut, rec_cut))

        # Standardized pAUC values (0.5..1)
        # For pAUROC: fpr_max is taken as the actual cutoff point for FPR,
        # or use the passed err_threshold for area consistency.
        e = float(fpr_cut[-1])
        pauroc = standardized_pauc(pauroc_raw, pauc_min=(e * e / 2.0), pauc_max=e)

        r = float(rec_cut[-1])
        pauprc = standardized_pauc(pauprc_raw, pauc_min=(0.5 * r), pauc_max=r)  # balanced baseline=0.5

        motif.statistics = {
                                "auPRC": auprc,
                                "auROC": auroc,
                                "pauPRC": pauprc,
                                "pauROC": pauroc,
                            }

        return {
            "PRC": {"RECALL": rec.tolist(), "PRECISION": prec.tolist()},
            "ROC": {"FPR": fpr.tolist(), "TPR": tpr.tolist()},
            "auPRC": auprc,
            "auROC": auroc,
            "pauPRC": pauprc,
            "pauROC": pauroc,
        }


class Bootstrapper:
    """Perform odd/even bootstrapping of motif discovery and evaluation.

    The bootstrapper splits the foreground sequences into two folds
    (odd and even indices), runs a motif discovery tool on the
    training fold and evaluates each discovered motif on the held out
    fold against a common background. The resulting statistics
    provide an unbiased estimate of motif generalisation performance.
    """

    def __init__(self, discovery_tool: MotifDiscoveryTool, evaluator: PerformanceEvaluator, output_dir: str) -> None:
        self.discovery_tool = discovery_tool
        self.evaluator = evaluator
        self.output_dir = output_dir

    def run(self, peaks: RaggedData, background: RaggedData,
            number_of_motifs: int, err_threshold: float, discovery_params: Dict[str, Iterable[Any]]
            ) -> Tuple[Dict[str, Any], List['MotifModel']]:
        """Execute bootstrap validation for motif discovery using odd/even sequence splitting.

        This method performs cross-validation by splitting the input sequences into odd and even
        indexed folds. For each fold, motif discovery is performed on the training set, and the
        discovered motifs are evaluated on the held-out test set. This process provides an
        unbiased estimate of motif generalization performance. Multiple parameter combinations
        are tested to optimize motif discovery settings.

        Parameters
        ----------
        peaks : RaggedData
            Foreground sequences to be split into training and test sets.
        background : RaggedData
            Background sequences used for motif discovery and evaluation.
        number_of_motifs : int
            Number of motifs to discover per bootstrap iteration.
        err_threshold : float
            FPR threshold for partial AUC computation during evaluation.
        discovery_params : Dict[str, Iterable[Any]]
            Dictionary mapping parameter names to lists of values to test.
            Example: {'length': [8, 10], 'order': [0, 1]}

        Returns
        -------
        statistics : Dict[str, Any]
            Dictionary mapping motif identifiers to their performance metrics.
        bootstrap_motifs : List[MotifModel]
            List of all discovered motif models across all bootstrap iterations.
        """
        statistics: Dict[str, Any] = {}
        bootstrap_motifs: List['MotifModel'] = []

        print("Starting bootstrap...")

        param_keys = list(discovery_params.keys())
        param_values = list(discovery_params.values())

        for combination in itertools.product(*param_values):
            current_params = dict(zip(param_keys, combination))
            params_suffix = format_params(current_params)

            for step_name in ["odd", "even"]:
                prefix_str = f"bootstrap_{params_suffix}_{step_name}"
                with tempfile.TemporaryDirectory(dir=os.path.join(self.output_dir, self.discovery_tool.name),
                                               prefix=prefix_str,
                                               delete=True) as tmp_dir:

                    # Split peaks into training and test based on index parity
                    n_peaks = peaks.num_sequences
                    if step_name == "odd":
                        train_indices = [i for i in range(n_peaks) if (i + 1) % 2 != 0]
                        test_indices = [i for i in range(n_peaks) if (i + 1) % 2 == 0]
                    else:
                        train_indices = [i for i in range(n_peaks) if (i + 1) % 2 == 0]
                        test_indices = [i for i in range(n_peaks) if (i + 1) % 2 != 0]

                    # For discovery tool we still need a list of arrays for write_fasta
                    # OR update write_fasta to handle RaggedData (which it already does).
                    # But discovery_tool.discover takes paths, so we must write files.
                    
                    train_peaks_list = [peaks.get_slice(i) for i in train_indices]
                    
                    # Create RaggedData for test set evaluation
                    test_peaks = ragged_from_list([peaks.get_slice(i) for i in test_indices], dtype=np.int8)

                    # Run motif discovery on the training peaks
                    fg_path = os.path.join(tmp_dir, "train.fasta")
                    bg_path = os.path.join(tmp_dir, "background.fasta")

                    # write_fasta accepts List[np.ndarray] OR RaggedData
                    write_fasta(train_peaks_list, fg_path)
                    write_fasta(background, bg_path)

                    motifs = self.discovery_tool.discover(
                        fg_path,
                        bg_path,
                        tmp_dir,
                        number_of_motifs=number_of_motifs,
                        **current_params
                    )

                    bootstrap_motifs += motifs

                    # Evaluate each motif on the test set
                    for motif in motifs:
                        stats = self.evaluator.evaluate(motif, test_peaks, background, err_threshold)
                        motif.name = f"{motif.name}_{params_suffix}_{step_name}"

                        statistics[motif.name] = stats

        return statistics, bootstrap_motifs
