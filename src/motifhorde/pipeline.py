"""Bootstrap, comparison, selection, and portable output orchestration."""

from __future__ import annotations

import itertools
import json
import os
import tempfile
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any, TypeVar

import numpy as np
from mimosa import (
    BestPerSequence,
    EncodedSequences,
    MotifModel,
    NullDistribution,
    TopFractionHits,
    read_fasta,
    reconstruct_pfm,
    write_model,
)

from .comparison import (
    ComparisonSettings,
    MatchedPair,
    comparison_value,
    default_threshold,
    deduplicate_pairs,
    filter_similar,
    is_similar,
    load_null_distribution,
    sort_pairs,
    validate_null_distribution,
)
from .discovery import MotifDiscoveryTool
from .evaluation import (
    EvaluatedMotif,
    Bootstrapper,
    PerformanceEvaluator,
    VALIDATION_METRICS,
)
from .io import write_meme

T = TypeVar("T")


@dataclass(frozen=True, slots=True)
class FinalMotif:
    """A full-data model promoted from a bootstrap pair with its metrics."""

    model: MotifModel
    name: str
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


class DeNovoPipeline:
    """Run the MotifHORDE parameter-selection workflow."""

    def __init__(
        self,
        discovery_tool: MotifDiscoveryTool,
        evaluator: PerformanceEvaluator,
        mimosa_metric: str = "co",
        comparison_criterion: str = "score",
        comparison_threshold: float | None = None,
        null_distribution: NullDistribution | str | os.PathLike | None = None,
        fpr_threshold: float = 0.001,
        number_of_motifs: int = 5,
        jobs: int = 1,
        seed: int | None = None,
        mimosa_search_range: int = 10,
    ) -> None:
        if number_of_motifs < 1:
            raise ValueError("number_of_motifs must be positive")
        if not np.isfinite(fpr_threshold) or not 0 <= fpr_threshold <= 1:
            raise ValueError("fpr_threshold must be finite and lie in [0, 1]")
        null = (
            load_null_distribution(null_distribution)
            if isinstance(null_distribution, (str, os.PathLike))
            else null_distribution
        )
        self.comparison = ComparisonSettings(
            metric=mimosa_metric,
            criterion=comparison_criterion,
            threshold=(
                default_threshold(comparison_criterion)
                if comparison_threshold is None
                else float(comparison_threshold)
            ),
            null_distribution=null,
            search_range=mimosa_search_range,
            jobs=jobs,
        )
        if comparison_criterion == "p-value":
            validate_null_distribution(
                null,
                metric=mimosa_metric,
                search_range=mimosa_search_range,
            )
        self.discovery_tool = discovery_tool
        self.evaluator = evaluator
        self.fpr_threshold = fpr_threshold
        self.number_of_motifs = number_of_motifs
        self.seed = seed

    def run(
        self,
        foreground_path: str,
        background_path: str,
        promoters_path: str,
        output_dir: str,
        discovery_params: dict[str, Iterable[Any]],
        metric: str,
    ) -> None:
        if metric not in VALIDATION_METRICS:
            raise ValueError(f"unsupported validation metric: {metric}")
        paths = self._prepare_output_dirs(output_dir)
        peaks, peak_names = self._read_sequences(foreground_path)
        background, background_names = self._read_sequences(background_path)
        promoters, _ = self._read_sequences(promoters_path)

        records = Bootstrapper(
            self.discovery_tool,
            self.evaluator,
            output_dir,
            jobs=self.comparison.jobs,
            seed=self.seed,
        ).run(
            peaks,
            background,
            self.number_of_motifs,
            self.fpr_threshold,
            discovery_params,
            peak_names,
            background_names,
        )
        self._save_bootstrap(records, paths["bootstrap"])
        match_groups = self._compare_bootstrap_motifs(
            records, discovery_params, peaks, metric
        )
        final_motifs = self._select_final_motifs(
            match_groups,
            peaks,
            promoters,
            foreground_path,
            background_path,
            output_dir,
            metric,
        )
        final_motifs = _select_nonredundant(
            [(motif.model, motif.metrics[metric], motif) for motif in final_motifs],
            peaks,
            self.comparison,
        )
        self._save_results(final_motifs, paths["motifs"], promoters, metric)

    def _prepare_output_dirs(self, output_dir: str) -> dict[str, Path]:
        root = Path(output_dir) / self.discovery_tool.name
        result = {
            "root": root,
            "bootstrap": root / "bootstrap",
            "motifs": root / "motifs",
        }
        conflicts = [
            result["bootstrap"] / "statistics.json",
            result["motifs"] / "statistics.json",
            result["motifs"] / "models" / "all_motifs_in_pfm_form.meme",
        ]
        for model_dir in (result["bootstrap"] / "models", result["motifs"] / "models"):
            if model_dir.exists() and any(model_dir.iterdir()):
                conflicts.append(model_dir)
        task_dir = result["bootstrap"] / "tasks"
        if task_dir.exists() and any(task_dir.iterdir()):
            conflicts.append(task_dir)
        existing = next((path for path in conflicts if path.exists()), None)
        if existing is not None:
            raise FileExistsError(f"MotifHORDE output already exists: {existing}")
        for path in result.values():
            path.mkdir(parents=True, exist_ok=True)
        return result

    @staticmethod
    def _read_sequences(path: str) -> tuple[EncodedSequences, tuple[str, ...]]:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"FASTA file not found: {path}")
        return read_fasta(path)

    def _save_bootstrap(self, records: list[EvaluatedMotif], result_dir: Path) -> None:
        models_dir = result_dir / "models"
        models_dir.mkdir(parents=True, exist_ok=True)
        entries = []
        paths = []
        for record in records:
            relative = Path("models") / _safe_name(record.id)
            target = result_dir / relative
            _check_new_path(target)
            paths.append(target)
            entries.append(
                {
                    "id": record.id,
                    "model": relative.as_posix(),
                    "name": record.model.name,
                    "split": record.split,
                    "params": dict(record.params),
                    "metrics": dict(record.metrics),
                }
            )
        for target, record in zip(paths, records, strict=True):
            write_model(target, record.model)
        _write_statistics(
            result_dir / "statistics.json", self.discovery_tool.name, entries
        )

    def _compare_bootstrap_motifs(
        self,
        records: list[EvaluatedMotif],
        discovery_params: dict[str, Iterable[Any]],
        peaks: EncodedSequences,
        metric: str,
    ) -> list[dict[str, Any]]:
        groups: list[dict[str, Any]] = []
        for params in _iter_param_grid(discovery_params):
            splits: list[list[EvaluatedMotif]] = []
            for split in ("odd", "even"):
                split_records = [
                    record
                    for record in records
                    if record.split == split and dict(record.params) == params
                ]
                if not split_records:
                    continue
                splits.append(
                    _select_nonredundant(
                        [
                            (record.model, record.metrics.get(metric, 0.0), record)
                            for record in split_records
                        ],
                        peaks,
                        self.comparison,
                    )
                )
            if len(splits) < 2:
                continue
            odd, even = splits
            pairs: list[MatchedPair] = []
            for query in odd:
                pairs.extend(self.comparison.compare(query, even, peaks))
            pairs = deduplicate_pairs(
                sort_pairs(
                    filter_similar(
                        pairs,
                        self.comparison.criterion,
                        self.comparison.threshold,
                    ),
                    self.comparison.criterion,
                )
            )
            if pairs:
                groups.append({"params": params, "pairs": pairs})
        return groups

    def _select_final_motifs(
        self,
        match_groups: list[dict[str, Any]],
        peaks: EncodedSequences,
        promoters: EncodedSequences,
        foreground_path: str,
        background_path: str,
        output_dir: str,
        metric: str,
    ) -> list[FinalMotif]:
        criterion = self.comparison.criterion
        threshold = self.comparison.threshold
        final: list[FinalMotif] = []
        for group in match_groups:
            params = dict(group["params"])
            with _temporary_discovery_dir(
                output_dir, self.discovery_tool.name
            ) as work_dir:
                full_motifs = self.discovery_tool.discover(
                    foreground_path,
                    background_path,
                    work_dir,
                    number_of_motifs=self.number_of_motifs * 2,
                    **params,
                )
            assigned: set[int] = set()
            pairs = sorted(
                group["pairs"],
                key=lambda pair: _average_metric(pair, metric),
                reverse=True,
            )
            for pair in pairs:
                candidates = []
                for index, model in enumerate(full_motifs):
                    if index in assigned:
                        continue
                    odd = self.comparison.compare(model, [pair.query], peaks)[0]
                    even = self.comparison.compare(model, [pair.target], peaks)[0]
                    if is_similar(odd, criterion, threshold) and is_similar(
                        even, criterion, threshold
                    ):
                        candidates.append(
                            (
                                (
                                    comparison_value(odd, criterion)
                                    + comparison_value(even, criterion)
                                )
                                / 2.0,
                                index,
                                model,
                            )
                        )
                if not candidates:
                    continue
                best = (
                    min(candidates, key=lambda item: item[0])
                    if criterion == "p-value"
                    else max(candidates, key=lambda item: item[0])
                )
                assigned.add(best[1])
                final.append(
                    FinalMotif(
                        model=best[2],
                        name=best[2].name,
                        params=params,
                        metrics=_average_metrics(pair),
                    )
                )
        return final

    def _save_results(
        self,
        records: list[FinalMotif],
        result_dir: Path,
        promoters: EncodedSequences,
        metric: str,
    ) -> None:
        ordered = sorted(records, key=lambda item: item.metrics[metric], reverse=True)
        models_dir = result_dir / "models"
        models_dir.mkdir(parents=True, exist_ok=True)
        paths = [
            models_dir / f"rank_{rank:03d}_{_safe_name(record.name)}"
            for rank, record in enumerate(ordered, start=1)
        ]
        meme_path = models_dir / "all_motifs_in_pfm_form.meme"
        stats_path = result_dir / "statistics.json"
        for path in [*paths, meme_path, stats_path]:
            _check_new_path(path)

        pfms = [_ensure_pfm(record.model, promoters) for record in ordered]
        write_meme(
            pfms,
            [(record.name, record.model.motif_length) for record in ordered],
            meme_path,
        )
        entries = []
        for rank, (path, record) in enumerate(
            zip(paths, ordered, strict=True), start=1
        ):
            write_model(path, record.model)
            entries.append(
                {
                    "rank": rank,
                    "model": path.relative_to(result_dir).as_posix(),
                    "name": record.name,
                    "params": dict(record.params),
                    "metrics": dict(record.metrics),
                }
            )
        _write_statistics(stats_path, self.discovery_tool.name, entries)


def _temporary_discovery_dir(output_dir: str, tool: str):
    return tempfile.TemporaryDirectory(
        prefix="full_", dir=os.path.join(output_dir, tool)
    )


def _check_new_path(path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")


def _write_statistics(path: Path, tool: str, motifs: list[dict[str, Any]]) -> None:
    path.write_text(
        json.dumps(
            {"schema_version": 1, "tool": tool, "motifs": motifs},
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )


def _ensure_pfm(model: MotifModel, sequences: EncodedSequences):
    return reconstruct_pfm(
        model,
        sequences,
        TopFractionHits(0.10, BestPerSequence()),
        strands="best",
    )


def _safe_name(name: str) -> str:
    return name.replace("/", "_").replace("\\", "_").replace(":", "-")


def _iter_param_grid(discovery_params: dict[str, Iterable[Any]]):
    keys = list(discovery_params)
    values = [tuple(discovery_params[key]) for key in keys]
    for combination in itertools.product(*values):
        yield dict(zip(keys, combination, strict=True))


def _average_metrics(pair: MatchedPair) -> dict[str, float]:
    return {
        metric: (float(pair.query.metrics[metric]) + float(pair.target.metrics[metric]))
        / 2.0
        for metric in VALIDATION_METRICS
    }


def _average_metric(pair: MatchedPair, metric: str) -> float:
    return (
        float(pair.query.metrics[metric]) + float(pair.target.metrics[metric])
    ) / 2.0


def _select_nonredundant(
    items: Sequence[tuple[MotifModel, float, T]],
    sequences: EncodedSequences,
    settings: ComparisonSettings,
) -> list[T]:
    """Greedily keep best-scoring items dissimilar to every kept item."""
    selected: list[tuple[MotifModel, float, T]] = []
    for item in sorted(items, key=lambda entry: entry[1], reverse=True):
        if selected and filter_similar(
            settings.compare(item[0], [kept[0] for kept in selected], sequences),
            settings.criterion,
            settings.threshold,
        ):
            continue
        selected.append(item)
    return [item[2] for item in selected]
