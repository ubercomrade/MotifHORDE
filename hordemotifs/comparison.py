"""
comparison
==========

Implementations of motif comparison metrics.  Comparing motifs is
useful for identifying similar patterns discovered in different
datasets or cross‑validation folds.  This module defines a common
interface for comparison algorithms and several concrete
implementations
"""

from __future__ import annotations

import os
import tempfile
from abc import ABC, abstractmethod
from typing import List, Optional


import numpy as np
import pandas as pd
from scipy.ndimage import convolve1d
from joblib import Parallel, delayed


from .models import MotifModel, RaggedScores
from .execute import run_motali_pwm_pwm
from .functions import (
    _fast_cj_kernel_numba,
    _fast_overlap_kernel_numba,
    _fast_pearson_kernel,
    scores_to_frequencies
)
from .ragged import RaggedData, ragged_from_list


class GeneralMotifComparator(ABC):
    """
    Abstract base class for motif comparators.
    
    This class defines the common interface for all motif comparison algorithms.
    Concrete implementations should inherit from this class and implement
    the compare method.
    """

    def __init__(self, name: str) -> None:
        """
        Initialize the comparator.
        
        Parameters
        ----------
        name : str
            Name of the comparator instance.
        """
        self.name = name

    @abstractmethod
    def compare(self,
                motifs_1: List[MotifModel],
                motifs_2: List[MotifModel],
                sequences: RaggedData | None = None,
                ) -> pd.DataFrame:
        """
        Compare motifs from two collections.
        
        This is an abstract method that must be implemented by subclasses.
        
        Parameters
        ----------
        motifs_1 : List[MotifModel]
            First collection of motifs to compare.
        motifs_2 : List[MotifModel]
            Second collection of motifs to compare.
        sequences : RaggedData or None
            Sequences for frequency calculation (if needed by the implementation).
            
        Returns
        -------
        pd.DataFrame
            DataFrame containing comparison results.
        """
        raise NotImplementedError


class TomtomComparator(GeneralMotifComparator):
    """
    Comparator for motifs using Euclidean Distance (ED) or Pearson Correlation (PCC).
    Includes Monte Carlo p-value estimation.
    """

    def __init__(
        self,
        metric: str = 'pcc',
        n_permutations: int = 1000,
        permute_rows: bool = False,
        pfm_mode: bool = False,
        n_jobs: int = 1,
        seed: Optional[int] = None
    ):
        """
        Initialize comparator.
        
        Parameters
        ----------
        metric : str
            'pcc' or 'ed'.
        n_permutations : int
            Number of Monte Carlo permutations for p-value calculation.
            Set to 0 to disable.
        permute_rows : bool
            If True, shuffles values within each column (destroys nucleotide structure).
            If False, only shuffles columns (positions).
        n_jobs : int
            Number of parallel jobs. -1 to use all cores.
        seed : int, optional
            Random seed for reproducibility.
        """
        super().__init__(name=f"TomtomComparator_{metric.upper()}")
        self.metric = metric.lower()
        if self.metric not in ['pcc', 'ed']:
            raise ValueError(f"Unsupported metric: {metric}. Use 'pcc' or 'ed'.")
        
        self.n_permutations = n_permutations
        self.permute_rows = permute_rows
        self.pfm_mode = pfm_mode
        self.n_jobs = n_jobs
        self.seed = seed

    def _prepare_matrix(self, matrix: np.ndarray):
        """
        Remove 'N' (if present), create Reverse Complement, and flatten nucleotide dimensions.
        
        Ensures consistent (4, L) shape for the new vectorized logic.
        
        Parameters
        ----------
        matrix : np.ndarray
            Input matrix in various formats (may include 'N' row).
            
        Returns
        -------
        tuple
            Tuple containing (original_matrix, reverse_complement_matrix) both in (4, L) format.
        """
        # 1. Remove 'N' row if present (assumes it's the last row)
        if matrix.shape[0] == 5:  # 5th row is for 'N'
            matrix = matrix[:-1]  # Keep only first 4 rows (A, C, G, T)
        
        # 2. Ensure matrix is in (4, L) format
        if matrix.shape[0] != 4:
            # If matrix is in (L, 4) format, transpose it
            if matrix.shape[1] == 4:
                matrix = matrix.T
            else:
                raise ValueError(f"Matrix has unexpected shape: {matrix.shape}. Expected (4, L) or (L, 4)")
        
        # 3. Create reverse complement
        # For PWM matrices: flip positions (axis=1) and bases (axis=0) to get RC
        # A->T, T->A, C->G, G->C: [A,C,G,T] becomes [T,G,C,A] when both axes flipped
        rc_matrix = np.flip(matrix, axis=0)  # Flip A/C/G/T
        rc_matrix = np.flip(rc_matrix, axis=1)  # Flip positions
        
        # 4. Return both original and reverse complement matrices
        # Both should be in (4, L) format
        return matrix, rc_matrix

    def _randomize_matrix(self, matrix: np.ndarray, rng: np.random.Generator):
        """
        Shuffle columns and optionally rows (values) in the original multidimensional matrix.
        
        This function implements a surrogate generation procedure where the nucleotide
        structure can be partially or completely destroyed depending on the permute_rows setting.
        
        Parameters
        ----------
        matrix : np.ndarray
            Input matrix to randomize.
        rng : np.random.Generator
            Random number generator instance.
            
        Returns
        -------
        np.ndarray
            Randomized matrix with shuffled columns and optionally rows.
        """
        # Work with a copy of the full dimensionality
        shuffled = matrix.copy()
        
        # 1. Shuffle columns (positions) along the last axis
        # Indices for the last axis
        pos_indices = np.arange(shuffled.shape[-1])
        rng.shuffle(pos_indices)
        shuffled = shuffled[..., pos_indices]
        
        # 2. Optional shuffle of "rows" (values within columns)
        # This destroys correlations (e.g., AA -> TG)
        if self.permute_rows:
            # Iterate through all positions
            for i in range(shuffled.shape[-1]):
                # Get slice for position i (e.g., 4x4 for dinucleotides)
                col_slice = shuffled[..., i]
                # Flatten, shuffle, reshape back
                flat_vals = col_slice.ravel()
                rng.shuffle(flat_vals)
                shuffled[..., i] = flat_vals.reshape(col_slice.shape)
                
        return shuffled

    def _vectorized_pcc(self, M1: np.ndarray, M2: np.ndarray):
        """
        Compute vectorized Pearson Correlation Coefficient between columns of M1 and M2.
        
        Parameters
        ----------
        M1 : np.ndarray
            Shape (4, L1) matrix representing first motif
        M2 : np.ndarray
            Shape (4, L2) matrix representing second motif
            
        Returns
        -------
        correlations : np.ndarray
            Array of correlations between corresponding columns
        """
        # Center both matrices by subtracting column means
        M1_centered = M1 - np.mean(M1, axis=0, keepdims=True)
        M2_centered = M2 - np.mean(M2, axis=0, keepdims=True)
        
        # Compute standard deviations for normalization
        M1_stds = np.sqrt(np.sum(M1_centered**2, axis=0))
        M2_stds = np.sqrt(np.sum(M2_centered**2, axis=0))
        
        # Handle zero-variance columns by setting std to 1 (will result in 0 correlation)
        M1_stds = np.where(M1_stds == 0, 1, M1_stds)
        M2_stds = np.where(M2_stds == 0, 1, M2_stds)
        
        # Compute dot product between centered matrices
        numerator = np.sum(M1_centered * M2_centered, axis=0)
        
        # Compute correlations
        denominators = M1_stds * M2_stds
        correlations = np.where(denominators != 0, numerator / denominators, 0.0)
        
        return correlations

    def _column_similarity(self, col1: np.ndarray, col2: np.ndarray):
        """
        Calculate similarity between two columns based on the selected metric.
        
        Parameters
        ----------
        col1 : np.ndarray
            First column for comparison.
        col2 : np.ndarray
            Second column for comparison.
            
        Returns
        -------
        float
            Similarity score based on the selected metric.
        """
        if self.metric == 'pcc':
            if np.std(col1) == 0 or np.std(col2) == 0:
                return 0.0
            return np.corrcoef(col1, col2)[0, 1]
        elif self.metric == 'ed':
            # Changed from Frobenius norm to sum of column-wise Euclidean distances
            dist = np.linalg.norm(col1 - col2)
            return -dist  # Negative because higher distance means less similarity
        return 0.0

    def _align_motifs(self, M1: np.ndarray, M2: np.ndarray):
        """
        Align two motifs by sliding one along the other and computing the best score.
        
        Parameters
        ----------
        M1 : np.ndarray
            First motif matrix of shape (4, L1).
        M2 : np.ndarray
            Second motif matrix of shape (4, L2).
            
        Returns
        -------
        tuple
            Tuple containing (best_score, best_offset) where:
            best_score : Best alignment score found.
            best_offset : Offset at which best score occurs.
        """
        L1 = M1.shape[1]
        L2 = M2.shape[1]
        best_score = -np.inf if self.metric == 'ed' else -np.inf
        best_offset = 0
        
        min_offset = -(L2 - 1)
        max_offset = L1 - 1
        
        for offset in range(min_offset, max_offset + 1):
            if offset < 0:
                len_overlap = min(L1, L2 + offset)
                if len_overlap <= 0:
                    continue
                cols1 = M1[:, :len_overlap]
                cols2 = M2[:, -offset : -offset + len_overlap]
            else:
                len_overlap = min(L1 - offset, L2)
                if len_overlap <= 0:
                    continue
                cols1 = M1[:, offset : offset + len_overlap]
                cols2 = M2[:, :len_overlap]
            
            if self.metric == 'ed':
                # Compute sum of column-wise Euclidean distances
                # This is the sum of ||col1_i - col2_i|| for each column pair
                column_distances = np.sqrt(np.sum((cols1 - cols2)**2, axis=0))
                current_score = -np.sum(column_distances)  # Negative because higher distance means less similarity
            elif self.metric == 'pcc':
                # Use vectorized PCC computation
                correlations = self._vectorized_pcc(cols1, cols2)
                current_score = np.sum(correlations)
            
            if current_score > best_score:
                best_score = current_score
                best_offset = offset
                
        return best_score, best_offset

    def _run_single_permutation(
        self,
        M1_flat: np.ndarray,
        M2_orig_matrix: np.ndarray,
        seed: int
    ):
        """
        Worker function for parallel execution.
        
        Generates one surrogate for M2 and compares it with M1.
        
        Parameters
        ----------
        M1_flat : np.ndarray
            Flattened version of the first motif matrix.
        M2_orig_matrix : np.ndarray
            Original matrix for the second motif (before flattening).
        seed : int
            Random seed for this permutation.
            
        Returns
        -------
        float
            Maximum alignment score between M1 and the randomized M2.
        """
        rng = np.random.default_rng(seed)
        
        # 1. Randomize the original M2 matrix (full dimensionality)
        M2_rand_matrix = self._randomize_matrix(M2_orig_matrix, rng)
        
        # 2. Prepare randomized matrix (flatten + rc)
        M2_rand_flat, M2_rand_rc_flat = self._prepare_matrix(M2_rand_matrix)
        
        # 3. Compare
        score_pp, _ = self._align_motifs(M1_flat, M2_rand_flat)
        score_pm, _ = self._align_motifs(M1_flat, M2_rand_rc_flat)
        
        return max(score_pp, score_pm)

    def compare(
        self,
        motifs_1: List[MotifModel],
        motifs_2: List[MotifModel],
        sequences: RaggedData | None = None,
    ) -> pd.DataFrame:
        """
        Compare lists of motifs with optional p-value calculation.
        
        Parameters
        ----------
        motifs_1 : List[MotifModel]
            List of first motif models to compare.
        motifs_2 : List[MotifModel]
            List of second motif models to compare.
            
        Returns
        -------
        pd.DataFrame
            DataFrame containing comparison results with columns:
            - query: name of the first motif
            - target: name of the second motif
            - score: alignment score
            - offset: optimal offset for alignment
            - orientation: strand orientation ('++', '+-')
            - metric: comparison metric used
            - p-value: statistical significance (if calculated)
            - z-score: standardized score (if calculated)
            - null_mean: mean of null distribution (if calculated)
            - null_std: std of null distribution (if calculated)
        """
        records = []
        if self.pfm_mode:
            for motif in {*motifs_1, *motifs_2}:
                motif.get_pfm(sequences=sequences, top_fraction=0.10)

        for m1 in motifs_1:
            # Prepare Query (M1)
            if self.pfm_mode:
                m1_flat, _ = self._prepare_matrix(m1.pfm)
            else:
                m1_flat, _ = self._prepare_matrix(m1.matrix)
            
            for m2 in motifs_2:
                # Prepare Target (M2)
                if self.pfm_mode:
                    m2_flat, m2_rc_flat = self._prepare_matrix(m2.pfm)
                else:
                    m2_flat, m2_rc_flat = self._prepare_matrix(m2.matrix)
                
                # --- Observed Score ---
                obs_score_pp, obs_off_pp = self._align_motifs(m1_flat, m2_flat)
                obs_score_pm, obs_off_pm = self._align_motifs(m1_flat, m2_rc_flat)
                
                if obs_score_pm > obs_score_pp:
                    obs_score = obs_score_pm
                    obs_offset = obs_off_pm
                    orientation = "+-"
                else:
                    obs_score = obs_score_pp
                    obs_offset = obs_off_pp
                    orientation = "++"
                
                row = {
                    "query": m1.name,
                    "target": m2.name,
                    "score": float(obs_score),
                    "offset": int(obs_offset),
                    "orientation": orientation,
                    "metric": self.metric
                }
                
                # --- Monte Carlo Permutations ---
                if self.n_permutations > 0:
                    base_rng = np.random.default_rng(self.seed)
                    seeds = base_rng.integers(0, 2**31, size=self.n_permutations)
                    
                    # Run in parallel
                    null_scores = Parallel(n_jobs=self.n_jobs, backend="loky")(
                        delayed(self._run_single_permutation)(
                            m1_flat, m2.matrix, int(seeds[i])
                        )
                        for i in range(self.n_permutations)
                    )
                    
                    null_scores = np.array(null_scores)
                    
                    # --- P-value calculation ---
                    # p = (count(null >= obs) + 1) / (N + 1)
                    # Fixed: ensure we count correctly where null_scores are greater than or equal to observed
                    n_ge = np.sum(null_scores >= obs_score)
                    p_value = (n_ge + 1.0) / (self.n_permutations + 1.0)
                    
                    # Additional statistics
                    mean_null = np.mean(null_scores)
                    std_null = np.std(null_scores)
                    z_score = (obs_score - mean_null) / (std_null + 1e-9)
                    
                    row.update({
                        "p-value": float(p_value),
                        "z-score": float(z_score),
                        "null_mean": float(mean_null),
                        "null_std": float(std_null)
                    })
                
                records.append(row)
                
        return pd.DataFrame(records)
    

class MotaliComparator(GeneralMotifComparator):
    """Comparator that wraps the Motali program.

    This comparator uses an external Motali program to compute similarity
    between Position Frequency Matrices (PFMs).
    """

    def __init__(self, fasta_path: str, threshold: float = 0.95, tmp_directory: str = '/tmp') -> None:
        """
        Initialize the MotaliComparator.
        
        Parameters
        ----------
        fasta_path : str
            Path to the FASTA file containing sequences for comparison.
        threshold : float, optional
            Minimum score threshold for filtering results (default is 0.95).
        tmp_directory : str, optional
            Directory for temporary files (default is '/tmp').
        """
        super().__init__(name="motali")
        self.threshold = threshold
        self.tmp_directory = tmp_directory
        self.fasta_path = fasta_path

    def compare(self, motifs_1: List[MotifModel], motifs_2: List[MotifModel], sequences: RaggedData | None = None) -> pd.DataFrame:
        """
        Compare motifs from two collections using the Motali program.
        
        Parameters
        ----------
        motifs_1 : List[MotifModel]
            First collection of motifs to compare.
        motifs_2 : List[MotifModel]
            Second collection of motifs to compare.
        sequences : RaggedData or None
            Sequences for comparison (not used in this implementation).
            
        Returns
        -------
        pd.DataFrame
            DataFrame containing comparison results with columns:
            - query: name of the first motif
            - target: name of the second motif
            - score: similarity score computed by Motali
        """
        records = list()
        for m1 in motifs_1:
            for m2 in motifs_2:
                with tempfile.TemporaryDirectory(dir=self.tmp_directory, delete=True) as tmp:
                    m1_path = os.path.join(tmp, "motif_1.pfm")
                    m2_path = os.path.join(tmp, "motif_2.pfm")

                    d1_path = os.path.join(tmp, "thresholds_1.dist")
                    d2_path = os.path.join(tmp, "thresholds_2.dist")

                    overlap_path = os.path.join(tmp, "overlap.txt")
                    all_path = os.path.join(tmp, "all.txt")
                    sta_path = os.path.join(tmp, "sta.txt")

                    m1.write_pfm(m1_path)
                    m2.write_pfm(m2_path)


                    m1.write_dist(d1_path)
                    m2.write_dist(d2_path)

                    score = run_motali_pwm_pwm(self.fasta_path,
                                               m1_path, m2_path,
                                               d1_path, d2_path,
                                               overlap_path,
                                               all_path,
                                               sta_path)
                    records.append({"query": m1.name, "target": m2.name, "score": score})
        records = pd.DataFrame(records)
        records = records[records['score'] >= self.threshold].reset_index(drop=True)

        return records


class UniversalMotifComparator(GeneralMotifComparator):
    """
    Universal comparator implementation that integrates functionality from both
    MotifComparator and CorrelationComparator. Supports Jaccard ('cj'),
    Overlap ('co'), and Pearson Correlation ('corr') metrics with optional
    permutation-based statistics.
    """

    def __init__(
        self,
        name: str = "UnifiedComparator",
        metric: str = 'cj',
        n_permutations: int = 1000,
        distortion_level: float = 0.4,
        n_jobs: int = -1,
        seed: Optional[int] = None,
        filter_type: Optional[str] = None,  # 'score' or 'p-value' or None
        filter_threshold: float = 0.05,
        min_kernel_size: int = 3,
        max_kernel_size: int = 11,
        search_range: int = 10
    ) -> None:
        """
        Initialize the unified comparator.

        Parameters
        ----------
        name : str
            Name of the comparator instance.
        metric : str
            Similarity metric to use: 'cj' (Continuous Jaccard),
            'co' (Continuous Overlap), or 'corr' (Pearson Correlation).
        n_permutations : int
            Number of permutations for statistical significance testing.
        distortion_level : float
            Level of distortion for surrogate generation (used for 'cj' and 'co').
        n_jobs : int
            Number of parallel jobs for permutations.
        seed : int, optional
            Random seed for reproducibility.
        filter_type : str, optional
            Type of filtering to apply to results ('score' or 'p-value').
        filter_threshold : float
            Threshold value for filtering.
        search_range : int
            Range to search for optimal offset alignment.
        """
        super().__init__(name)
        self.metric = metric.lower()
        if self.metric not in ['cj', 'co', 'corr']:
            raise ValueError(f"Unsupported metric: {metric}. Use 'cj', 'co', or 'corr'.")

        self.n_permutations = n_permutations
        self.distortion_level = distortion_level
        self.n_jobs = n_jobs
        self.seed = seed
        self.filter_type = filter_type
        self.filter_threshold = filter_threshold
        self.min_kernel_size = min_kernel_size
        self.max_kernel_size = max_kernel_size
        self.search_range = search_range

    @staticmethod
    def _compute_metric_internal(S1: RaggedData, S2: RaggedData, search_range: int, metric: str):
        """Internal dispatcher for metric computation kernels."""
        if metric == "cj":
            return _fast_cj_kernel_numba(S1.data, S1.offsets, S2.data, S2.offsets, search_range)
        elif metric == "co":
            return _fast_overlap_kernel_numba(S1.data, S1.offsets, S2.data, S2.offsets, search_range)
        elif metric == "corr":
            # Returns (correlation, p-value, offset)
            return _fast_pearson_kernel(S1.data, S1.offsets, S2.data, S2.offsets, search_range)
        else:
            raise ValueError(f"Unknown metric: {metric}")

    def _single_compare(self, motif1: MotifModel, motif2: MotifModel, sequences: RaggedData, stat: bool = True):
        """Perform a single comparison between two motifs."""
        freq1_plus = motif1.get_frequencies(sequences, strand="+")
        freq2_plus = motif2.get_frequencies(sequences, strand="+")
        freq2_minus = motif2.get_frequencies(sequences, strand="-")

        # Observed scores for both orientations
        res_pp = self._compute_metric_internal(freq1_plus, freq2_plus, self.search_range, self.metric)
        res_pm = self._compute_metric_internal(freq1_plus, freq2_minus, self.search_range, self.metric)

        # Extract scores for comparison (first element of return tuple for all kernels)
        score_pp = res_pp[0]
        score_pm = res_pm[0]

        if score_pm > score_pp:
            orientation = "+-"
            obs_res = res_pm
            freq1, freq2 = freq1_plus, freq2_minus
        else:
            orientation = "++"
            obs_res = res_pp
            freq1, freq2 = freq1_plus, freq2_plus

        obs_score = float(obs_res[0])
        obs_offset = int(obs_res[-1])  # Offset is always the last element

        result = {
            "score": obs_score,
            "offset": obs_offset,
            "orientation": orientation,
            "strand_pair": orientation,
            "metric": self.metric
        }

        # Handle 'corr' specific p-value if stats are not requested via permutations
        if self.metric == 'corr':
            result["p-value"] = float(obs_res[1])

        if stat and self.n_permutations > 0:
            base_rng = np.random.default_rng(self.seed)
            seeds = base_rng.integers(0, 2**31, size=self.n_permutations)

            # Use MotifComparator's surrogate logic for 'cj' and 'co'
            # For 'corr', we use the same surrogate logic if permutations are requested
            results = Parallel(n_jobs=self.n_jobs, backend="loky")(
                delayed(self._compute_surrogate_score)(
                    freq1,
                    freq2,
                    np.random.default_rng(int(seeds[i])),
                )
                for i in range(self.n_permutations)
            )

            null_scores = np.array([r[0] for r in results], dtype=np.float32)
            null_mean = float(np.mean(null_scores))
            null_std = float(np.std(null_scores))
            n_ge = int(np.sum(null_scores >= obs_score))

            result.update({
                "p-value": (n_ge + 1.0) / (self.n_permutations + 1.0),
                "z-score": (obs_score - null_mean) / (null_std + 1e-9),
                "null_mean": null_mean,
                "null_std": null_std,
            })

        return result

    def _compute_surrogate_score(self, freq1: RaggedData, freq2: RaggedData, rng: np.random.Generator):
        """Helper for parallel permutation execution."""
        # Reusing the static method from MotifComparator as requested by "preserving algorithmic nuances"
        surrogate = self._generate_single_surrogate(
            freq2, 
            rng, 
            min_kernel_size=self.min_kernel_size, 
            max_kernel_size=self.max_kernel_size, 
            distortion_level=self.distortion_level
        )
        return self._compute_metric_internal(freq1, surrogate, self.search_range, self.metric)

    @staticmethod
    def _generate_single_surrogate(
        frequencies: RaggedData,
        rng: np.random.Generator,
        min_kernel_size: int = 3,
        max_kernel_size: int = 11,
        distortion_level: float = 1.0,
    ) -> RaggedData:
        """
        Generate a single surrogate frequency profile using convolution with a distorted kernel.
        
        This function implements a sophisticated surrogate generation algorithm that creates
        distorted versions of the input frequency profiles. The "distortion" logic refers to
        how the identity kernel is systematically modified through several techniques:
        
        1. Base kernel selection (smooth, edge, double_peak patterns)
        2. Noise addition with controlled amplitude
        3. Gradient application to introduce directional bias
        4. Smoothing to reduce artifacts
        5. Convex combination with identity kernel based on distortion level
        6. Sign flipping for additional variation
        
        Parameters
        ----------
        frequencies : RaggedData
            Input frequency profile to generate surrogate from.
        rng : np.random.Generator
            Random number generator instance.
        min_kernel_size : int, optional
            Minimum size of the convolution kernel (default is 3).
        max_kernel_size : int, optional
            Maximum size of the convolution kernel (default is 11).
        distortion_level : float, optional
            Level of distortion to apply (0.0 to 1.0, default is 1.0).
            
        Returns
        -------
        RaggedData
            Surrogate frequency profile generated from the input.
        """
        # For simplicity in surrogate generation, we use dense adapter
        dense_adapter = RaggedScores.from_numba(frequencies)
        X = dense_adapter.values
        lengths = dense_adapter.lengths

        kernel_size = int(rng.integers(min_kernel_size, max_kernel_size + 1))
        if kernel_size % 2 == 0:
            kernel_size += 1
        center = kernel_size // 2

        kernel_types = ["smooth", "edge", "double_peak"]
        kernel_type = str(rng.choice(kernel_types))

        identity_kernel = np.zeros(kernel_size, dtype=np.float32)
        identity_kernel[center] = 1.0

        if kernel_type == "smooth":
            x = np.linspace(-3, 3, kernel_size)
            base = np.exp(-0.5 * x**2).astype(np.float32)
        elif kernel_type == "edge":
            base = np.zeros(kernel_size, dtype=np.float32)
            base[max(center - 1, 0)] = -1.0
            base[min(center + 1, kernel_size - 1)] = 1.0
        elif kernel_type == "double_peak":
            base = np.zeros(kernel_size, dtype=np.float32)
            base[0] = 0.5
            base[-1] = 0.5
            base[center] = -1.0
        else:
            base = identity_kernel.copy()

        noise = rng.normal(0, 1, size=kernel_size).astype(np.float32)
        slope = float(rng.uniform(-1.0, 1.0)) * distortion_level * 2.0
        gradient = np.linspace(-slope, slope, kernel_size).astype(np.float32)

        distorted_kernel = base + distortion_level * noise + gradient

        if kernel_size >= 3:
            smooth_filter = np.array([0.25, 0.5, 0.25], dtype=np.float32)
            distorted_kernel = np.convolve(distorted_kernel, smooth_filter, mode="same")

        distorted_kernel /= (np.linalg.norm(distorted_kernel) + 1e-8)

        alpha = max(0.0, min(1.0, distortion_level))
        final_kernel = (1.0 - alpha) * identity_kernel + alpha * distorted_kernel
        if rng.uniform() < 0.5:
            final_kernel = -final_kernel
        final_kernel /= (np.linalg.norm(final_kernel) + 1e-8)

        convolved = convolve1d(X, final_kernel, axis=1, mode="constant", cval=0.0).astype(np.float32)
        
        # Convert back to RaggedData
        convolved_list = [convolved[i, :lengths[i]] for i in range(len(lengths))]
        convolved_ragged = ragged_from_list(convolved_list, dtype=np.float32)

        return scores_to_frequencies(convolved_ragged)

    def compare(
        self,
        motifs_1: List[MotifModel],
        motifs_2: List[MotifModel],
        sequences: RaggedData | None = None
    ) -> pd.DataFrame:
        """Compare motifs from two collections using the unified logic."""
        if sequences is None:
            raise ValueError("Sequences (RaggedData) are required for this comparator.")

        records = []
        calc_stats = self.n_permutations > 0

        for m1 in motifs_1:
            for m2 in motifs_2:
                out = self._single_compare(m1, m2, sequences, stat=calc_stats)
                row = {"query": m1.name, "target": m2.name}
                row.update(out)
                records.append(row)

        df = pd.DataFrame(records)

        if self.filter_type == 'score':
            df = df[df['score'] >= self.filter_threshold].reset_index(drop=True)
        elif self.filter_type == 'p-value' and 'p-value' in df.columns:
            df = df[df['p-value'] <= self.filter_threshold].reset_index(drop=True)

        return df
