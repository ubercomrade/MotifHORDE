"""Bootstrap, comparison, selection, and portable output orchestration."""

from __future__ import annotations

import itertools
import json
import os
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import numpy as np
from mimosa import (
    BestPerSequence,
    EncodedSequences,
    MotifModel,
    TopFractionHits,
    read_fasta,
    reconstruct_pfm,
    write_model,
)

from .comparison import (
    MatchedPair,
    SUPPORTED_METRICS,
    compare_records,
    deduplicate_pairs,
    default_threshold,
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


class DeNovoPipeline:
    """Run the MotifHORDE parameter-selection workflow."""

    def __init__(
        self,
        discovery_tool: MotifDiscoveryTool,
        evaluator: PerformanceEvaluator,
        mimosa_metric: str = "co",
        comparison_criterion: str = "score",
        comparison_threshold: float | None = None,
        null_distribution=None,
        fpr_threshold: float = 0.001,
        number_of_motifs: int = 5,
        jobs: int = 1,
        seed: int | None = None,
        mimosa_search_range: int = 10,
    ) -> None:
        if comparison_criterion not in {"score", "p-value"}:
            raise ValueError("comparison criterion must be 'score' or 'p-value'")
        if mimosa_metric not in SUPPORTED_METRICS:
            raise ValueError(f"unsupported Mimosa metric: {mimosa_metric}")
        if jobs < 1:
            raise ValueError("jobs must be positive")
        if mimosa_search_range < 0:
            raise ValueError("mimosa_search_range must be non-negative")
        if number_of_motifs < 1:
            raise ValueError("number_of_motifs must be positive")
        if not np.isfinite(fpr_threshold) or not 0 <= fpr_threshold <= 1:
            raise ValueError("fpr_threshold must be finite and lie in [0, 1]")
        if comparison_threshold is not None and not np.isfinite(comparison_threshold):
            raise ValueError("comparison threshold must be finite")
        self.discovery_tool = discovery_tool
        self.evaluator = evaluator
        self.mimosa_metric = mimosa_metric
        self.comparison_criterion = comparison_criterion
        self.comparison_threshold = (
            default_threshold(comparison_criterion)
            if comparison_threshold is None
            else float(comparison_threshold)
        )
        self.null_distribution = (
            load_null_distribution(null_distribution)
            if isinstance(null_distribution, (str, os.PathLike))
            else null_distribution
        )
        if comparison_criterion == "p-value" and self.null_distribution is None:
            raise ValueError("p-value comparison requires a Mimosa null bundle")
        if comparison_criterion == "p-value":
            validate_null_distribution(
                self.null_distribution,
                metric=mimosa_metric,
                search_range=mimosa_search_range,
            )
        self.fpr_threshold = fpr_threshold
        self.number_of_motifs = number_of_motifs
        self.jobs = jobs
        self.seed = seed
        self.mimosa_search_range = mimosa_search_range

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
            jobs=self.jobs,
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
        final_records = self._select_final_motifs(
            match_groups,
            peaks,
            promoters,
            foreground_path,
            background_path,
            output_dir,
            metric,
        )
        final_records = _deduplicate_final_motifs(
            final_records,
            metric,
            peaks,
            self.mimosa_metric,
            self.jobs,
            self.comparison_criterion,
            self.comparison_threshold,
            self.null_distribution,
            self.mimosa_search_range,
        )
        self._save_results(final_records, paths["motifs"], promoters, metric)

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
            odd = [
                record
                for record in records
                if record.split == "odd" and dict(record.params) == params
            ]
            even = [
                record
                for record in records
                if record.split == "even" and dict(record.params) == params
            ]
            if not odd or not even:
                continue
            odd = _select_nonredundant_motifs(
                odd,
                metric,
                peaks,
                self.mimosa_metric,
                self.jobs,
                self.comparison_criterion,
                self.comparison_threshold,
                self.null_distribution,
                self.mimosa_search_range,
            )
            even = _select_nonredundant_motifs(
                even,
                metric,
                peaks,
                self.mimosa_metric,
                self.jobs,
                self.comparison_criterion,
                self.comparison_threshold,
                self.null_distribution,
                self.mimosa_search_range,
            )
            pairs: list[MatchedPair] = []
            for query in odd:
                pairs.extend(
                    compare_records(
                        query,
                        even,
                        peaks,
                        metric=self.mimosa_metric,
                        jobs=self.jobs,
                        null_distribution=self.null_distribution,
                        search_range=self.mimosa_search_range,
                    )
                )
            pairs = deduplicate_pairs(
                sort_pairs(
                    filter_similar(
                        pairs,
                        self.comparison_criterion,
                        self.comparison_threshold,
                    ),
                    self.comparison_criterion,
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
    ) -> list[dict[str, Any]]:
        final: list[dict[str, Any]] = []
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
                    odd = compare_records(
                        model,
                        [pair.query],
                        peaks,
                        metric=self.mimosa_metric,
                        jobs=self.jobs,
                        null_distribution=self.null_distribution,
                        search_range=self.mimosa_search_range,
                    )[0]
                    even = compare_records(
                        model,
                        [pair.target],
                        peaks,
                        metric=self.mimosa_metric,
                        jobs=self.jobs,
                        null_distribution=self.null_distribution,
                        search_range=self.mimosa_search_range,
                    )[0]
                    if is_similar(
                        odd, self.comparison_criterion, self.comparison_threshold
                    ) and is_similar(
                        even, self.comparison_criterion, self.comparison_threshold
                    ):
                        candidates.append(
                            (
                                _candidate_value(odd, even, self.comparison_criterion),
                                index,
                                model,
                            )
                        )
                if not candidates:
                    continue
                best = (
                    min(candidates, key=lambda item: item[0])
                    if self.comparison_criterion == "p-value"
                    else max(candidates, key=lambda item: item[0])
                )
                assigned.add(best[1])
                final.append(
                    {
                        "model": best[2],
                        "name": best[2].name,
                        "params": params,
                        "metrics": _average_metrics(pair),
                    }
                )
        return final

    def _save_results(
        self,
        records: list[dict[str, Any]],
        result_dir: Path,
        promoters: EncodedSequences,
        metric: str,
    ) -> None:
        ordered = sorted(
            records, key=lambda item: item["metrics"][metric], reverse=True
        )
        models_dir = result_dir / "models"
        models_dir.mkdir(parents=True, exist_ok=True)
        paths = [
            models_dir / f"rank_{rank:03d}_{_safe_name(record['name'])}"
            for rank, record in enumerate(ordered, start=1)
        ]
        meme_path = models_dir / "all_motifs_in_pfm_form.meme"
        stats_path = result_dir / "statistics.json"
        for path in [*paths, meme_path, stats_path]:
            _check_new_path(path)

        pfms = [_ensure_pfm(record["model"], promoters) for record in ordered]
        write_meme(
            pfms,
            [(record["name"], record["model"].motif_length) for record in ordered],
            meme_path,
        )
        entries = []
        for rank, (path, record) in enumerate(
            zip(paths, ordered, strict=True), start=1
        ):
            write_model(path, record["model"])
            entries.append(
                {
                    "rank": rank,
                    "model": path.relative_to(result_dir).as_posix(),
                    "name": record["name"],
                    "params": dict(record["params"]),
                    "metrics": dict(record["metrics"]),
                }
            )
        _write_statistics(stats_path, self.discovery_tool.name, entries)


def _temporary_discovery_dir(output_dir: str, tool: str):
    import tempfile

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
    return _average_metrics(pair)[metric]


def _select_nonredundant_motifs(
    motifs: list[EvaluatedMotif],
    metric: str,
    sequences: EncodedSequences,
    mimosa_metric: str,
    jobs: int,
    criterion: str,
    threshold: float,
    null_distribution,
    search_range: int,
) -> list[EvaluatedMotif]:
    selected: list[EvaluatedMotif] = []
    for motif in sorted(
        motifs, key=lambda item: item.metrics.get(metric, 0.0), reverse=True
    ):
        if not selected:
            selected.append(motif)
            continue
        pairs = compare_records(
            motif,
            selected,
            sequences,
            metric=mimosa_metric,
            jobs=jobs,
            null_distribution=null_distribution,
            search_range=search_range,
        )
        if not filter_similar(pairs, criterion, threshold):
            selected.append(motif)
    return selected


def _candidate_value(first: MatchedPair, second: MatchedPair, criterion: str) -> float:
    from .comparison import comparison_value

    return (
        comparison_value(first, criterion) + comparison_value(second, criterion)
    ) / 2.0


def _deduplicate_final_motifs(
    records: list[dict[str, Any]],
    metric: str,
    sequences: EncodedSequences,
    mimosa_metric: str,
    jobs: int,
    criterion: str,
    threshold: float,
    null_distribution,
    search_range: int,
) -> list[dict[str, Any]]:
    kept: list[dict[str, Any]] = []
    for record in sorted(
        records, key=lambda item: item["metrics"][metric], reverse=True
    ):
        if not kept:
            kept.append(record)
            continue
        pairs = compare_records(
            record["model"],
            [item["model"] for item in kept],
            sequences,
            metric=mimosa_metric,
            jobs=jobs,
            null_distribution=null_distribution,
            search_range=search_range,
        )
        if not filter_similar(pairs, criterion, threshold):
            kept.append(record)
    return kept
