"""De novo motif discovery pipeline.

High level orchestration of the motif discovery workflow. The
pipeline performs an odd/even bootstrap to estimate motif quality,
compares motifs discovered in each fold, combines statistics and
selects the best motifs for further de‑novo discovery on the full
dataset. This implementation builds on the modular discovery,
comparison and evaluation classes defined elsewhere in the package.

"""
from __future__ import annotations

import os
import json
import tempfile
import itertools
from typing import Dict, List, Tuple, Any, Iterable, Optional

import pandas as pd
from .io import read_fasta, write_meme
from .models import MotifModel
from .discovery import MotifDiscoveryTool
from .evaluation import PerformanceEvaluator, Bootstrapper
from .comparison import UniversalMotifComparator, TomtomComparator, MotaliComparator
from .functions import format_params
from .ragged import RaggedData


class DeNovoPipeline:
    """De novo motif discovery pipeline combining bootstrapping, motif comparison and final selection.

    This pipeline implements a robust de novo motif discovery approach that uses bootstrap
    resampling to assess motif reliability and select high-confidence motifs for further analysis.
    
    Parameters
    ----------
    discovery_tool : MotifDiscoveryTool
        Tool used for motif discovery (e.g., STREME, MEME, DREME).
    evaluator : PerformanceEvaluator
        Tool used for evaluating motif performance.
    comparator : UniversalMotifComparator
        Tool used for comparing motifs (e.g., TOMTOM, custom comparator).
    fpr_threshold : float, optional
        False positive rate threshold for motif filtering (default is 0.001).
    number_of_motifs : int, optional
        Number of motifs to discover in each bootstrap iteration (default is 5).

    Attributes
    ----------
    discovery_tool : MotifDiscoveryTool
        Tool used for motif discovery.
    evaluator : PerformanceEvaluator
        Tool used for evaluating motif performance.
    comparator : UniversalMotifComparator
        Tool used for comparing motifs.
    fpr_threshold : float
        False positive rate threshold for motif filtering.
    number_of_motifs : int
        Number of motifs to discover in each bootstrap iteration.
    """

    def __init__(
        self,
        discovery_tool: MotifDiscoveryTool,
        evaluator: PerformanceEvaluator,
        comparator: UniversalMotifComparator,
        fpr_threshold: float = 0.001,
        number_of_motifs: int = 5,
    ) -> None:
        self.discovery_tool = discovery_tool
        self.evaluator = evaluator
        self.comparator = comparator
        self.fpr_threshold = fpr_threshold
        self.number_of_motifs = number_of_motifs

    def run(
        self,
        foreground_path: str,
        background_path: str,
        promoters_path: str,
        output_dir: str,
        discovery_params: Dict[str, Iterable[Any]],
        metric: str
    ) -> None:
        """Execute the complete de novo motif discovery pipeline.

        The pipeline performs the following high-level steps:
        
        1. Bootstrap motif discovery: Perform odd/even bootstrap sampling to identify
           reproducible motifs.
        2. Motif comparison: Compare motifs discovered in odd vs even bootstrap samples.
        3. Full dataset discovery: Discover motifs on the complete dataset.
        4. Motif selection: Select the best motifs based on bootstrap stability and
           performance metrics.
        5. Results saving: Save selected motifs and associated statistics.

        Parameters
        ----------
        foreground_path : str
            Path to the foreground FASTA file containing sequences of interest.
        background_path : str
            Path to the background FASTA file for statistical comparison.
        promoters_path : str
            Path to the promoters FASTA file for motif evaluation.
        output_dir : str
            Directory to store pipeline outputs.
        discovery_params : Dict[str, Iterable[Any]]
            Dictionary mapping parameter names to their possible values for discovery.
        metric : str
            Name of the metric to use for motif ranking (e.g., 'auPRC', 'auROC').
            
        Returns
        -------
        None
            Outputs are written to the specified output directory.
        """
        # Prepare output directories for bootstrap and final motif results
        bootstrap_dir, motifs_dir = self._prepare_output_dirs(output_dir)

        # Read input sequences (foreground peaks, background sequences, and promoters)
        peaks, background, promoters = self._read_sequences(
            foreground_path, background_path, promoters_path
        )

        # 1. Execute bootstrap discovery to assess motif reproducibility
        statistics, bootstrap_motifs = self._run_bootstrap(
            peaks, background, discovery_params, output_dir
        )

        # 2. Save bootstrap motifs with safe filenames to prevent collisions
        bootstrap_models_dir = os.path.join(bootstrap_dir, "models")
        os.makedirs(bootstrap_models_dir, exist_ok=True)

        print(f"Saving {len(bootstrap_motifs)} bootstrap motifs to {bootstrap_models_dir}...")
        for i, motif in enumerate(bootstrap_motifs):
            safe_name = motif.name.replace("/", "_").replace("\\", "_").replace(":", "-")
            # Add index to prevent name collisions if names happen to match
            filename = f"{i:04d}_{safe_name}.pkl"
            motif.save(os.path.join(bootstrap_models_dir, filename), clear_cache=True)

        # Save raw statistics from bootstrap runs
        self._save_json(statistics, os.path.join(bootstrap_dir, "statistics.json"))

        # 3. Compute PFM matrices for bootstrap motifs and compare odd/even pairs
        for motif in bootstrap_motifs:
            if motif.pfm is None:
                motif.get_pfm(peaks)

        bootstrap_records = self._compare_bootstrap_motifs(
            bootstrap_motifs, discovery_params, peaks, promoters, statistics
        )

        if bootstrap_records is None:
            print("No motif comparisons were made; exiting.")
            return

        # 4. Select final motifs based on full dataset discovery and bootstrap validation
        final_motifs, final_info, final_stats = self._select_final_motifs(
            bootstrap_records=bootstrap_records,
            bootstrap_motifs=bootstrap_motifs,
            peaks=peaks,
            promoters=promoters,
            metric=metric,
            foreground_path=foreground_path,
            background_path=background_path,
            output_dir=output_dir,
            discovery_params=discovery_params  # Pass parameters for proper configuration
        )

        # 5. Save final results including motifs, metadata, and performance statistics
        self._save_results(
             final_motifs, final_info, final_stats, motifs_dir, metric
        )

    def _prepare_output_dirs(self, output_dir: str) -> Tuple[str, str]:
        """Create and organize output directories for pipeline results.

        Parameters
        ----------
        output_dir : str
            Base directory for all pipeline outputs.

        Returns
        -------
        Tuple[str, str]
            Paths to bootstrap and motifs subdirectories respectively.
        """
        os.makedirs(output_dir, exist_ok=True)
        bootstrap_dir = os.path.join(output_dir, self.discovery_tool.name, "bootstrap")
        motifs_dir = os.path.join(output_dir, self.discovery_tool.name, "motifs")
        os.makedirs(bootstrap_dir, exist_ok=True)
        os.makedirs(motifs_dir, exist_ok=True)

        return bootstrap_dir, motifs_dir

    def _read_sequences(
        self, foreground_path: str, background_path: str, promoters_path: str
    ) -> Tuple[RaggedData, RaggedData, RaggedData]:
        """Read input sequences from FASTA files.

        Parameters
        ----------
        foreground_path : str
            Path to the foreground FASTA file.
        background_path : str
            Path to the background FASTA file.
        promoters_path : str
            Path to the promoters FASTA file.

        Returns
        -------
        Tuple[RaggedData, RaggedData, RaggedData]
            RaggedData objects for peaks, background, and promoters respectively.
        """
        # Validate input file paths exist and are accessible
        if not os.path.exists(foreground_path):
            raise FileNotFoundError(f"Foreground file not found: {foreground_path}")
        if not os.path.exists(background_path):
            raise FileNotFoundError(f"Background file not found: {background_path}")
        if not os.path.exists(promoters_path):
            raise FileNotFoundError(f"Promoters file not found: {promoters_path}")
            
        peaks = read_fasta(foreground_path, return_ragged=True)
        background = read_fasta(background_path, return_ragged=True)
        promoters = read_fasta(promoters_path, return_ragged=True)

        # Validate that sequences were loaded
        if peaks.num_sequences == 0:
            raise ValueError(f"No sequences found in foreground file: {foreground_path}")
        if background.num_sequences == 0:
            raise ValueError(f"No sequences found in background file: {background_path}")
        if promoters.num_sequences == 0:
            raise ValueError(f"No sequences found in promoters file: {promoters_path}")

        return peaks, background, promoters

    def _run_bootstrap(
        self, peaks: RaggedData, background: RaggedData, discovery_params: Dict[str, Iterable[Any]], output_dir: str
    ):
        """Execute bootstrap resampling for motif discovery.

        Parameters
        ----------
        peaks : RaggedData
            Foreground sequences for motif discovery.
        background : RaggedData
            Background sequences for statistical comparison.
        discovery_params : Dict[str, Iterable[Any]]
            Parameter combinations to test during discovery.
        output_dir : str
            Output directory for bootstrap results.

        Returns
        -------
        Tuple[Dict, List[MotifModel]]
            Statistics dictionary and list of discovered bootstrap motifs.
        """
        bootstrapper = Bootstrapper(self.discovery_tool, self.evaluator, output_dir)

        return bootstrapper.run(
            peaks, background, self.number_of_motifs, self.fpr_threshold, discovery_params
        )

    def _compare_bootstrap_motifs(
        self,
        bootstrap_motifs: List[MotifModel],
        discovery_params: Dict[str, Iterable[Any]],
        peaks: RaggedData,
        promoters: RaggedData,
        statistics: Dict,
    ) -> Optional[pd.DataFrame]:
        """Compare bootstrap motifs and generate combined statistics.

        This method compares motifs discovered in odd/even bootstrap samples
        for each parameter combination and combines performance statistics.

        Parameters
        ----------
        bootstrap_motifs : List[MotifModel]
            List of motifs discovered during bootstrap iterations.
        discovery_params : Dict[str, Iterable[Any]]
            Parameter combinations used for discovery.
        peaks : RaggedData
            Peak sequences for motif comparison.
        promoters : RaggedData
            Promoter sequences for motif evaluation.
        statistics : Dict
            Performance statistics from bootstrap runs.

        Returns
        -------
        Optional[pd.DataFrame]
            DataFrame with motif comparison results, or None if no comparisons made.
        """
        bootstrap_records = []

        comparator_sequences = peaks
        if isinstance(self.comparator, MotaliComparator):
            for motif in bootstrap_motifs:
                motif.get_threshold_table(promoters)

        param_keys = sorted(discovery_params.keys())
        param_values_list = [discovery_params[k] for k in param_keys]

        for combination in itertools.product(*param_values_list):
            current_params = dict(zip(param_keys, combination))
            param_suffix = format_params(current_params)

            odd_suffix = f"_{param_suffix}_odd"
            even_suffix = f"_{param_suffix}_even"

            odd_motifs = [m for m in bootstrap_motifs if m.name.endswith(odd_suffix)]
            even_motifs = [m for m in bootstrap_motifs if m.name.endswith(even_suffix)]

            if not odd_motifs or not even_motifs:
                continue

            record = self.comparator.compare(odd_motifs, even_motifs, sequences=comparator_sequences)

            for key, value in current_params.items():
                record[key] = value

            bootstrap_records.append(record)

        if not bootstrap_records:
            return None

        df = pd.concat(bootstrap_records, ignore_index=True)

        if "p-value" in df.columns:
            df = df.sort_values(by="p-value", ascending=True)
        elif "score" in df.columns:
            df = df.sort_values(by="score", ascending=False)
        else:
            raise ValueError("Comparison must contain either 'p-value' or 'score'.")

        # Remove duplicate entries
        df = df.drop_duplicates(subset=['query'])
        df = df.drop_duplicates(subset=['target'])
        df = df.reset_index(drop=True)

        # Combine statistics by averaging across folds
        for metric in ["auPRC", "auROC", "pauPRC", "pauROC"]:
            df[metric] = [
                (statistics.get(row["query"], {}).get(metric, 0.0) +
                 statistics.get(row["target"], {}).get(metric, 0.0)) / 2.0
                for _, row in df.iterrows()
            ]
        return df

    def _select_final_motifs(
        self,
        bootstrap_records: pd.DataFrame,
        bootstrap_motifs: List[MotifModel],
        peaks: RaggedData,
        promoters: RaggedData,
        metric: str,
        foreground_path: str,
        background_path: str,
        output_dir: str,
        discovery_params: Dict[str, Iterable[Any]],
    ):
        """Select final motifs based on bootstrap validation and full dataset discovery.

        This method performs de novo discovery on the full dataset and matches
        the resulting motifs to the most stable motifs identified through bootstrap.

        Parameters
        ----------
        bootstrap_records : pd.DataFrame
            Comparison results from bootstrap motif analysis.
        bootstrap_motifs : List[MotifModel]
            Motifs discovered during bootstrap iterations.
        peaks : RaggedData
            Peak sequences for motif evaluation.
        promoters : RaggedData
            Promoter sequences for motif evaluation.
        metric : str
            Metric to use for motif ranking.
        foreground_path : str
            Path to foreground sequences.
        background_path : str
            Path to background sequences.
        output_dir : str
            Output directory for results.
        discovery_params : Dict[str, Iterable[Any]]
            Parameter combinations for discovery.

        Returns
        -------
        Tuple[List[MotifModel], List[Tuple[str, Dict[str, Any]]], Dict[str, Dict[str, float]]]
            Final motifs, metadata, and performance statistics.
        """
        final_motifs: List[MotifModel] = []
        final_info: List[Tuple[str, Dict[str, Any]]] = []
        final_stats: Dict[str, Dict[str, float]] = {}

        param_keys = sorted(discovery_params.keys())

        # Process groupby for single or multiple parameters
        for group_key, group in bootstrap_records.groupby(param_keys):
            # group_key can be: 12  or (12,) or (12, 2, ...)
            if not isinstance(group_key, tuple):
                group_key = (group_key,)

            # Handle rare edge case of tuple inside tuple
            group_key = tuple(v[0] if isinstance(v, tuple) and len(v) == 1 else v for v in group_key)

            current_params = dict(zip(param_keys, group_key))
            param_suffix = format_params(current_params)

            # Sort group by quality metric
            if metric in group.columns:
                group = group.sort_values(metric, ascending=False)
            else:
                raise ValueError(f"Comparison must contain {metric}")

            with tempfile.TemporaryDirectory(dir=os.path.join(output_dir, self.discovery_tool.name)) as tmp_dir:
                # Run discovery on full dataset
                full_motifs_list = self.discovery_tool.discover(
                    foreground_path,
                    background_path,
                    tmp_dir,
                    number_of_motifs=self.number_of_motifs * 2,
                    **current_params
                )
                for motif in full_motifs_list:
                    if motif.pfm is None:
                        motif.get_pfm(promoters)

                assigned = set()
                for _, record in group.iterrows():
                    odd_name, even_name = record["query"], record["target"]

                    remaining_fulls = [m for m in full_motifs_list if m.name not in assigned]
                    if not remaining_fulls:
                        break

                    best_full_motif = self._select_best_full_motif(
                        remaining_fulls, odd_name, even_name, bootstrap_motifs, peaks, promoters,
                    )

                    log_prefix = f"Params {current_params}:"

                    if best_full_motif is None:
                        print(f"{log_prefix} No match found for motifs {odd_name}, {even_name}")
                        continue
                    else:
                        print(f"{log_prefix} Best match for {odd_name} and {even_name} is {best_full_motif.name}")


                    assigned.add(best_full_motif.name)

                    # Save info
                    final_motifs.append(best_full_motif)
                    final_info.append((best_full_motif.name, current_params))
                    stat_key = f"{best_full_motif.name}_{param_suffix}"
                    final_stats[stat_key] = {
                        m: record[m] for m in ["auPRC", "auROC", "pauPRC", "pauROC"]
                    }

        return final_motifs, final_info, final_stats

    def _select_best_full_motif(
        self,
        full_motifs: List[MotifModel],
        odd_name: str,
        even_name: str,
        bootstrap_motifs: List[MotifModel],
        peaks: RaggedData,
        promoters: RaggedData,
    ) -> Optional[MotifModel]:
        """Match full dataset motifs to bootstrap-validated motifs.

        This method compares motifs discovered on the full dataset to the
        corresponding odd/even bootstrap motifs and selects the best matching motif.

        Parameters
        ----------
        full_motifs : List[MotifModel]
            Motifs discovered on the full dataset.
        odd_name : str
            Name of the odd bootstrap motif.
        even_name : str
            Name of the even bootstrap motif.
        bootstrap_motifs : List[MotifModel]
            Motifs discovered during bootstrap iterations.
        peaks : RaggedData
            Peak sequences for motif comparison.
        promoters : RaggedData
            Promoter sequences for motif evaluation.

        Returns
        -------
        Optional[MotifModel]
            Best matching motif from the full dataset discovery, or None if no match found.
        """
        # preparation
        comparator_sequences = None if isinstance(self.comparator, TomtomComparator) else peaks
        if isinstance(self.comparator, MotaliComparator):
            for motif in full_motifs:
                motif.get_threshold_table(promoters)
            for motif in bootstrap_motifs:
                motif.get_threshold_table(promoters)


        odd_ref = [m for m in bootstrap_motifs if m.name == odd_name]
        even_ref = [m for m in bootstrap_motifs if m.name == even_name]

        comparison_odd = self.comparator.compare(full_motifs, odd_ref, sequences=comparator_sequences)
        comparison_even = self.comparator.compare(full_motifs, even_ref, sequences=comparator_sequences)

        if "p-value" in comparison_odd.columns:
            metric, selector = "p-value", min
        elif "score" in comparison_odd.columns:
            metric, selector = "score", max
        else:
            raise ValueError("Comparison must contain either 'p-value' or 'score'.")

        avg_scores = {}
        for motif in full_motifs:
            # Filter by motif name (query)
            odd_vals = comparison_odd.loc[comparison_odd["query"] == motif.name, metric].values
            even_vals = comparison_even.loc[comparison_even["query"] == motif.name, metric].values

            if len(odd_vals) == 0 or len(even_vals) == 0:
                continue

            # Average the metric (score or p-value)
            avg_scores[motif.name] = (odd_vals[0] + even_vals[0]) / 2

        if not avg_scores:
            return None

        # Select the best by metric (min for p-value, max for score)
        best_name = selector(avg_scores, key=avg_scores.get)
        return next((m for m in full_motifs if m.name == best_name), None)

    def _save_results(
        self,
        final_motifs: List[MotifModel],
        final_info: List[Tuple[str, Dict[str, Any]]],
        final_stats: Dict[str, Dict[str, float]],
        motifs_dir: str,
        metric: str
    ) -> None:
        """Save final motifs, metadata, and performance statistics.

        This method organizes and saves the final results of the pipeline,
        including ranked motifs, metadata, and performance statistics.

        Parameters
        ----------
        final_motifs : List[MotifModel]
            Selected motifs to save.
        final_info : List[Tuple[str, Dict[str, Any]]]
            Metadata for each motif (name and parameters).
        final_stats : Dict[str, Dict[str, float]]
            Performance statistics for each motif.
        motifs_dir : str
            Directory to save motif results.
        metric : str
            Primary metric used for ranking.

        Returns
        -------
        None
            Results are written to disk.
        """
        def get_stats_key(name: str, params: Dict[str, Any]) -> str:
            return f"{name}_{format_params(params)}"

        # 1. Sort indices by selected metric
        sorted_indices = sorted(
            range(len(final_info)),
            key=lambda i: final_stats[get_stats_key(final_info[i][0], final_info[i][1])][metric],
            reverse=True,
        )

        motifs_sorted = [final_motifs[i] for i in sorted_indices]
        info_sorted = [final_info[i] for i in sorted_indices]

        # 2. Models directory
        models_output_dir = os.path.join(motifs_dir, "models")
        os.makedirs(models_output_dir, exist_ok=True)
        print(f"Saving {len(motifs_sorted)} individual models to {models_output_dir}...")

        # 3.1 MEME file with all PFMs
        pfms = [m.pfm for m in motifs_sorted]
        meta_data = [(m.name, m.length) for m in motifs_sorted]
        write_meme(pfms, meta_data, f"{models_output_dir}/all_motifs_in_pfm_form.meme")

        # 3.2 Save models
        for rank, motif in enumerate(motifs_sorted, start=1):
            safe_name = motif.name.replace("/", "_").replace("\\", "_").replace(":", "-")
            filename = f"{rank:03d}_{safe_name}.pkl"
            motif.save(os.path.join(models_output_dir, filename), clear_cache=True)

        # 4. statistics.json
        self._save_json(final_stats, os.path.join(motifs_dir, "statistics.json"))

        # 5. Print all metrics to console
        metrics_to_print = ["auPRC", "auROC", "pauPRC", "pauROC"]

        for idx, (name, params) in enumerate(info_sorted, start=1):
            key = get_stats_key(name, params)
            stats = final_stats[key]  # expected dict with metrics
            params_str = ", ".join(f"{k}={params[k]}" for k in sorted(params.keys()))

            # Print selected metric + others (if available)
            parts = [f"{m}={stats[m]:.4f}" for m in metrics_to_print if m in stats]
            print(f"Motif {idx}: {name}; {params_str}; " + "; ".join(parts))


    @staticmethod
    def _save_json(data: Dict, path: str) -> None:
        """Save data to JSON file with indentation.

        Parameters
        ----------
        data : Dict
            Data to serialize to JSON.
        path : str
            Path to the output JSON file.

        Returns
        -------
        None
            Data is written to the specified file.
        """
        with open(path, "w") as handle:
            json.dump(data, handle, indent=2)
