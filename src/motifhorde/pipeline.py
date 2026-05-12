"""De novo motif discovery pipeline."""

from __future__ import annotations

import itertools
import json
import os
import tempfile
from typing import Any, Dict, Iterable, List, Optional, Tuple

import joblib
import pandas as pd

from .batches import SequenceBatch
from .comparison import TomtomComparator, UniversalMotifComparator
from .discovery import MotifDiscoveryTool
from .evaluation import Bootstrapper, PerformanceEvaluator
from .functions import format_params
from .io import read_fasta, write_meme
from .models import GenericModel, get_pfm

SIMILARITY_PVALUE_THRESHOLD = 0.001
SIMILARITY_SCORE_THRESHOLD = 0.9
VALIDATION_METRICS = ("auPRC", "auROC", "pauPRC", "pauROC")


class DeNovoPipeline:
    """Orchestrate bootstrap validation, motif comparison, and final discovery."""

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
        metric: str,
    ) -> None:
        bootstrap_dir, motifs_dir = self._prepare_output_dirs(output_dir)
        peaks, background, promoters = self._read_sequences(foreground_path, background_path, promoters_path)

        statistics, bootstrap_motifs = self._run_bootstrap(peaks, background, discovery_params, output_dir)
        self._save_bootstrap(bootstrap_motifs, statistics, bootstrap_dir)

        bootstrap_records = self._compare_bootstrap_motifs(bootstrap_motifs, discovery_params, peaks, statistics, metric)
        if bootstrap_records is None:
            print("No motif comparisons were made; exiting.")
            return

        final_motifs, final_info, final_stats = self._select_final_motifs(
            bootstrap_records=bootstrap_records,
            bootstrap_motifs=bootstrap_motifs,
            peaks=peaks,
            promoters=promoters,
            metric=metric,
            foreground_path=foreground_path,
            background_path=background_path,
            output_dir=output_dir,
            discovery_params=discovery_params,
        )
        final_motifs, final_info, final_stats = _deduplicate_final_motifs(
            final_motifs,
            final_info,
            final_stats,
            metric,
            self.comparator,
            peaks,
        )
        self._save_results(final_motifs, final_info, final_stats, motifs_dir, metric, promoters)

    def _prepare_output_dirs(self, output_dir: str) -> Tuple[str, str]:
        os.makedirs(output_dir, exist_ok=True)
        bootstrap_dir = os.path.join(output_dir, self.discovery_tool.name, "bootstrap")
        motifs_dir = os.path.join(output_dir, self.discovery_tool.name, "motifs")
        os.makedirs(bootstrap_dir, exist_ok=True)
        os.makedirs(motifs_dir, exist_ok=True)
        return bootstrap_dir, motifs_dir

    def _read_sequences(
        self,
        foreground_path: str,
        background_path: str,
        promoters_path: str,
    ) -> Tuple[SequenceBatch, SequenceBatch, SequenceBatch]:
        for label, path in {
            "Foreground": foreground_path,
            "Background": background_path,
            "Promoters": promoters_path,
        }.items():
            if not os.path.exists(path):
                raise FileNotFoundError(f"{label} file not found: {path}")

        peaks = read_fasta(foreground_path)
        background = read_fasta(background_path)
        promoters = read_fasta(promoters_path)

        for label, batch, path in (
            ("foreground", peaks, foreground_path),
            ("background", background, background_path),
            ("promoters", promoters, promoters_path),
        ):
            if len(batch["lengths"]) == 0:
                raise ValueError(f"No sequences found in {label} file: {path}")

        return peaks, background, promoters

    def _run_bootstrap(
        self,
        peaks: SequenceBatch,
        background: SequenceBatch,
        discovery_params: Dict[str, Iterable[Any]],
        output_dir: str,
    ):
        bootstrapper = Bootstrapper(self.discovery_tool, self.evaluator, output_dir)
        return bootstrapper.run(peaks, background, self.number_of_motifs, self.fpr_threshold, discovery_params)

    def _save_bootstrap(self, motifs: List[GenericModel], statistics: Dict[str, Any], bootstrap_dir: str) -> None:
        models_dir = os.path.join(bootstrap_dir, "models")
        os.makedirs(models_dir, exist_ok=True)
        print(f"Saving {len(motifs)} bootstrap motifs to {models_dir}...")
        for index, motif in enumerate(motifs):
            joblib.dump(motif, os.path.join(models_dir, f"{index:04d}_{_safe_name(motif.name)}.pkl"))
        self._save_json(statistics, os.path.join(bootstrap_dir, "statistics.json"))

    def _compare_bootstrap_motifs(
        self,
        bootstrap_motifs: List[GenericModel],
        discovery_params: Dict[str, Iterable[Any]],
        peaks: SequenceBatch,
        statistics: Dict[str, Any],
        metric: str,
    ) -> Optional[pd.DataFrame]:
        records = []

        for current_params in _iter_param_grid(discovery_params):
            odd_motifs, even_motifs = _bootstrap_motifs_for_params(bootstrap_motifs, current_params)
            if not odd_motifs or not even_motifs:
                continue

            odd_selected = _select_nonredundant_motifs(odd_motifs, statistics, metric, self.comparator, peaks)
            even_selected = _select_nonredundant_motifs(even_motifs, statistics, metric, self.comparator, peaks)
            if not odd_selected or not even_selected:
                continue

            frame = self.comparator.compare(odd_selected, even_selected, sequences=peaks)
            frame = _filter_similar_matches(frame)
            if frame.empty:
                continue
            frame = _deduplicate_matches(_sort_comparisons(frame))
            for key, value in current_params.items():
                frame[key] = value
            frame = _attach_average_metrics(frame, statistics)
            records.append(frame)

        if not records:
            return None

        return pd.concat(records, ignore_index=True).reset_index(drop=True)

    def _select_final_motifs(
        self,
        bootstrap_records: pd.DataFrame,
        bootstrap_motifs: List[GenericModel],
        peaks: SequenceBatch,
        promoters: SequenceBatch,
        metric: str,
        foreground_path: str,
        background_path: str,
        output_dir: str,
        discovery_params: Dict[str, Iterable[Any]],
    ):
        final_motifs: List[GenericModel] = []
        final_info: List[Tuple[str, Dict[str, Any]]] = []
        final_stats: Dict[str, Dict[str, float]] = {}
        param_keys = sorted(discovery_params.keys())

        for group_key, group in bootstrap_records.groupby(param_keys):
            if not isinstance(group_key, tuple):
                group_key = (group_key,)
            current_params = dict(zip(param_keys, group_key))
            param_suffix = format_params(current_params)

            if metric not in group.columns:
                raise ValueError(f"Comparison must contain {metric}")
            group = group.sort_values(metric, ascending=False)

            with tempfile.TemporaryDirectory(dir=os.path.join(output_dir, self.discovery_tool.name)) as tmp_dir:
                full_motifs = self.discovery_tool.discover(
                    foreground_path,
                    background_path,
                    tmp_dir,
                    number_of_motifs=self.number_of_motifs * 2,
                    **current_params,
                )

                assigned = set()
                for _, record in group.iterrows():
                    remaining = [motif for motif in full_motifs if motif.name not in assigned]
                    if not remaining:
                        break

                    best = self._select_best_full_motif(
                        remaining,
                        record["query"],
                        record["target"],
                        bootstrap_motifs,
                        peaks,
                    )
                    if best is None:
                        print(f"Params {current_params}: No match found for motifs {record['query']}, {record['target']}")
                        continue

                    print(f"Params {current_params}: Best match for {record['query']} and {record['target']} is {best.name}")
                    assigned.add(best.name)
                    _ensure_pfm(best, promoters)
                    final_motifs.append(best)
                    final_info.append((best.name, current_params))
                    final_stats[f"{best.name}_{param_suffix}"] = {
                        name: record[name] for name in ["auPRC", "auROC", "pauPRC", "pauROC"]
                    }

        return final_motifs, final_info, final_stats

    def _select_best_full_motif(
        self,
        full_motifs: List[GenericModel],
        odd_name: str,
        even_name: str,
        bootstrap_motifs: List[GenericModel],
        peaks: SequenceBatch,
    ) -> Optional[GenericModel]:
        odd_ref = [motif for motif in bootstrap_motifs if motif.name == odd_name]
        even_ref = [motif for motif in bootstrap_motifs if motif.name == even_name]
        if not odd_ref or not even_ref:
            return None

        sequences = None if isinstance(self.comparator, TomtomComparator) else peaks
        comparison_odd = self.comparator.compare(full_motifs, odd_ref, sequences=sequences)
        comparison_even = self.comparator.compare(full_motifs, even_ref, sequences=sequences)

        compare_metric = _comparison_column(comparison_odd)
        _comparison_column(comparison_even)

        candidates = []
        for motif in full_motifs:
            odd_values = comparison_odd.loc[comparison_odd["query"] == motif.name, compare_metric].values
            even_values = comparison_even.loc[comparison_even["query"] == motif.name, compare_metric].values
            if not len(odd_values) or not len(even_values):
                continue
            odd_value = float(odd_values[0])
            even_value = float(even_values[0])
            if not _is_similar_value(compare_metric, odd_value):
                continue
            if not _is_similar_value(compare_metric, even_value):
                continue
            candidates.append((motif, (odd_value + even_value) / 2.0))

        if not candidates:
            return None
        return sorted(candidates, key=lambda item: item[1], reverse=not _comparison_sort_ascending(compare_metric))[0][0]

    def _save_results(
        self,
        final_motifs: List[GenericModel],
        final_info: List[Tuple[str, Dict[str, Any]]],
        final_stats: Dict[str, Dict[str, float]],
        motifs_dir: str,
        metric: str,
        promoters: SequenceBatch,
    ) -> None:
        def stats_key(name: str, params: Dict[str, Any]) -> str:
            return f"{name}_{format_params(params)}"

        sorted_indices = sorted(
            range(len(final_info)),
            key=lambda index: final_stats[stats_key(final_info[index][0], final_info[index][1])][metric],
            reverse=True,
        )
        motifs_sorted = [final_motifs[index] for index in sorted_indices]
        info_sorted = [final_info[index] for index in sorted_indices]

        models_dir = os.path.join(motifs_dir, "models")
        os.makedirs(models_dir, exist_ok=True)
        print(f"Saving {len(motifs_sorted)} individual models to {models_dir}...")

        pfms = [_ensure_pfm(motif, promoters) for motif in motifs_sorted]
        metadata = [(motif.name, motif.length) for motif in motifs_sorted]
        write_meme(pfms, metadata, os.path.join(models_dir, "all_motifs_in_pfm_form.meme"))

        for rank, motif in enumerate(motifs_sorted, start=1):
            joblib.dump(motif, os.path.join(models_dir, f"{rank:03d}_{_safe_name(motif.name)}.pkl"))

        self._save_json(final_stats, os.path.join(motifs_dir, "statistics.json"))

        for index, (name, params) in enumerate(info_sorted, start=1):
            stats = final_stats[stats_key(name, params)]
            params_str = ", ".join(f"{key}={params[key]}" for key in sorted(params))
            parts = [f"{key}={stats[key]:.4f}" for key in ["auPRC", "auROC", "pauPRC", "pauROC"] if key in stats]
            print(f"Motif {index}: {name}; {params_str}; " + "; ".join(parts))

    @staticmethod
    def _save_json(data: Dict, path: str) -> None:
        with open(path, "w") as handle:
            json.dump(data, handle, indent=2)


def _ensure_pfm(model: GenericModel, sequences: SequenceBatch):
    pfm = model.config.get("_source_pfm")
    if pfm is None:
        pfm = model.config.get("_derived_pfm")
    if pfm is None:
        pfm = get_pfm(model, sequences, top_fraction=0.10)
        model.config["_derived_pfm"] = pfm
    return pfm


def _safe_name(name: str) -> str:
    return name.replace("/", "_").replace("\\", "_").replace(":", "-")


def _iter_param_grid(discovery_params: Dict[str, Iterable[Any]]) -> Iterable[Dict[str, Any]]:
    param_keys = sorted(discovery_params)
    param_values = [discovery_params[key] for key in param_keys]
    for combination in itertools.product(*param_values):
        yield dict(zip(param_keys, combination))


def _bootstrap_motifs_for_params(
    bootstrap_motifs: List[GenericModel],
    current_params: Dict[str, Any],
) -> Tuple[List[GenericModel], List[GenericModel]]:
    param_suffix = format_params(current_params)
    odd_suffix = f"_{param_suffix}_odd"
    even_suffix = f"_{param_suffix}_even"
    odd_motifs = [motif for motif in bootstrap_motifs if motif.name.endswith(odd_suffix)]
    even_motifs = [motif for motif in bootstrap_motifs if motif.name.endswith(even_suffix)]
    return odd_motifs, even_motifs


def _comparison_column(frame: pd.DataFrame) -> str:
    if "p-value" in frame.columns:
        return "p-value"
    if "score" in frame.columns:
        return "score"
    raise ValueError("Comparison must contain either 'p-value' or 'score'.")


def _comparison_sort_ascending(column: str) -> bool:
    if column == "p-value":
        return True
    if column == "score":
        return False
    raise ValueError("Comparison column must be either 'p-value' or 'score'.")


def _is_similar_value(column: str, value: float) -> bool:
    if pd.isna(value):
        return False
    if column == "p-value":
        return value <= SIMILARITY_PVALUE_THRESHOLD
    if column == "score":
        return value >= SIMILARITY_SCORE_THRESHOLD
    raise ValueError("Comparison column must be either 'p-value' or 'score'.")


def _filter_similar_matches(frame: pd.DataFrame) -> pd.DataFrame:
    column = _comparison_column(frame)
    mask = frame[column].map(lambda value: _is_similar_value(column, float(value)))
    return frame[mask].reset_index(drop=True)


def _sort_comparisons(frame: pd.DataFrame) -> pd.DataFrame:
    column = _comparison_column(frame)
    return frame.sort_values(by=column, ascending=_comparison_sort_ascending(column)).reset_index(drop=True)


def _deduplicate_matches(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.drop_duplicates(subset=["query"]).drop_duplicates(subset=["target"]).reset_index(drop=True)


def _select_nonredundant_motifs(
    motifs: List[GenericModel],
    statistics: Dict[str, Any],
    metric: str,
    comparator,
    sequences,
) -> List[GenericModel]:
    sorted_motifs = sorted(
        motifs,
        key=lambda motif: statistics.get(motif.name, {}).get(metric, 0.0),
        reverse=True,
    )
    selected: List[GenericModel] = []
    for motif in sorted_motifs:
        if not selected:
            selected.append(motif)
            continue

        frame = comparator.compare([motif], selected, sequences=sequences)
        similar = _filter_similar_matches(frame)
        if similar.empty:
            selected.append(motif)
    return selected


def _attach_average_metrics(frame: pd.DataFrame, statistics: Dict[str, Any]) -> pd.DataFrame:
    frame = frame.copy()
    for metric in VALIDATION_METRICS:
        frame[metric] = [
            (statistics.get(row["query"], {}).get(metric, 0.0) + statistics.get(row["target"], {}).get(metric, 0.0))
            / 2.0
            for _, row in frame.iterrows()
        ]
    return frame


def _deduplicate_final_motifs(
    final_motifs: List[GenericModel],
    final_info: List[Tuple[str, Dict[str, Any]]],
    final_stats: Dict[str, Dict[str, float]],
    metric: str,
    comparator,
    sequences,
) -> Tuple[List[GenericModel], List[Tuple[str, Dict[str, Any]]], Dict[str, Dict[str, float]]]:
    def stats_key(name: str, params: Dict[str, Any]) -> str:
        return f"{name}_{format_params(params)}"

    sorted_indices = sorted(
        range(len(final_motifs)),
        key=lambda index: final_stats[stats_key(final_info[index][0], final_info[index][1])][metric],
        reverse=True,
    )

    kept_indices: List[int] = []
    kept_motifs: List[GenericModel] = []
    for index in sorted_indices:
        motif = final_motifs[index]
        if not kept_motifs:
            kept_indices.append(index)
            kept_motifs.append(motif)
            continue

        frame = comparator.compare([motif], kept_motifs, sequences=sequences)
        similar = _filter_similar_matches(frame)
        if similar.empty:
            kept_indices.append(index)
            kept_motifs.append(motif)

    deduplicated_motifs = [final_motifs[index] for index in kept_indices]
    deduplicated_info = [final_info[index] for index in kept_indices]
    deduplicated_stats = {
        stats_key(name, params): final_stats[stats_key(name, params)] for name, params in deduplicated_info
    }
    return deduplicated_motifs, deduplicated_info, deduplicated_stats
