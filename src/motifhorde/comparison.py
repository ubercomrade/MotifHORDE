"""Pipeline-facing comparison wrappers backed by MIMOSA."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

import pandas as pd
from mimosa.api import compare_one_to_many as mimosa_compare_one_to_many
from mimosa.comparison import (
    SUPPORTED_MOTIF_METRICS,
    SUPPORTED_PROFILE_METRICS,
    compare,
    compare_one_to_many as low_level_compare_one_to_many,
    create_comparator_config,
)
from mimosa.types import ComparatorConfig, ComparisonResult

DEFAULT_SCORE_COMPARISON_THRESHOLD = 0.9
DEFAULT_PVALUE_COMPARISON_THRESHOLD = 0.05
MISSING_ADJUSTED_PVALUE_MESSAGE = (
    "Comparison results do not contain adjusted p-values. Provide a compatible "
    "MIMOSA null distribution with --mimosa-null-distribution."
)


def comparison_column_for_criterion(criterion: str) -> str:
    if criterion == "score":
        return "score"
    if criterion == "p-value":
        return "adj.p-value"
    raise ValueError(f"Unsupported comparison criterion: {criterion}")


def default_threshold_for_criterion(criterion: str) -> float:
    if criterion == "score":
        return DEFAULT_SCORE_COMPARISON_THRESHOLD
    if criterion == "p-value":
        return DEFAULT_PVALUE_COMPARISON_THRESHOLD
    raise ValueError(f"Unsupported comparison criterion: {criterion}")


class GeneralMotifComparator:
    """Small pipeline-facing comparator interface."""

    def __init__(
        self,
        name: str,
        comparison_criterion: str = "score",
        comparison_threshold: float | None = None,
    ) -> None:
        self.name = name
        self.comparison_criterion = comparison_criterion
        self.comparison_threshold = (
            default_threshold_for_criterion(comparison_criterion)
            if comparison_threshold is None
            else comparison_threshold
        )
        self.comparison_column = comparison_column_for_criterion(comparison_criterion)

    def compare(self, motifs_1, motifs_2, sequences=None) -> pd.DataFrame:
        raise NotImplementedError


class TomtomComparator(GeneralMotifComparator):
    """Matrix comparator compatible with the pipeline naming."""

    def __init__(
        self,
        metric: str = "pcc",
        pfm_mode: bool = False,
        n_jobs: int | None = 1,
        seed: int | None = None,
        comparison_criterion: str = "score",
        comparison_threshold: float | None = None,
        *,
        pvalue: bool = False,
        null_distribution: str | Path | dict[str, Any] | None = None,
        null_search_dirs: Iterable[str | Path] | None = None,
        effective_number_of_targets: int | None = None,
    ) -> None:
        super().__init__(
            name=f"TomtomComparator_{metric.upper()}",
            comparison_criterion=comparison_criterion,
            comparison_threshold=comparison_threshold,
        )
        self.config = create_comparator_config(
            metric=metric,
            pfm_mode=pfm_mode,
            n_jobs=n_jobs,
            seed=seed,
            pvalue=pvalue,
            null_distribution=null_distribution,
            null_search_dirs=null_search_dirs,
            effective_number_of_targets=effective_number_of_targets,
        )

    def compare(self, motifs_1, motifs_2, sequences=None) -> pd.DataFrame:
        records = []
        for motif in motifs_1:
            records.extend(
                low_level_compare_one_to_many(
                    motif, motifs_2, "motif", self.config, sequences=sequences
                )
            )
        return _records_to_frame(records)


class UniversalMotifComparator(GeneralMotifComparator):
    """Profile comparator wrapper for the de novo pipeline."""

    def __init__(
        self,
        name: str = "UnifiedComparator",
        metric: str = "co",
        n_jobs: int | None = -1,
        seed: int | None = None,
        comparison_criterion: str = "score",
        comparison_threshold: float | None = None,
        search_range: int = 10,
        min_logfpr: float | None = None,
        window_radius: int = 10,
        realign_window: int = 3,
        cache_mode: str = "off",
        cache_dir: str = ".mimosa-cache",
        *,
        pvalue: bool = False,
        null_distribution: str | Path | dict[str, Any] | None = None,
        null_search_dirs: Iterable[str | Path] | None = None,
        effective_number_of_targets: int | None = None,
    ) -> None:
        super().__init__(
            name,
            comparison_criterion=comparison_criterion,
            comparison_threshold=comparison_threshold,
        )
        self.config = create_comparator_config(
            metric=metric,
            n_jobs=n_jobs,
            seed=seed,
            search_range=search_range,
            min_logfpr=min_logfpr,
            window_radius=window_radius,
            realign_window=realign_window,
            cache_mode=cache_mode,
            cache_dir=cache_dir,
            pvalue=pvalue,
            null_distribution=null_distribution,
            null_search_dirs=null_search_dirs,
            effective_number_of_targets=effective_number_of_targets,
        )

    def compare(self, motifs_1, motifs_2, sequences=None) -> pd.DataFrame:
        records = []
        for motif in motifs_1:
            records.extend(
                mimosa_compare_one_to_many(
                    query=motif,
                    targets=list(motifs_2),
                    strategy="profile",
                    sequences=sequences,
                    comparator=self.config,
                )
            )
        return _records_to_frame(records)


def _records_to_frame(records):
    return pd.DataFrame.from_records(records)


compare_one_to_many = low_level_compare_one_to_many


__all__ = [
    "ComparatorConfig",
    "ComparisonResult",
    "DEFAULT_PVALUE_COMPARISON_THRESHOLD",
    "DEFAULT_SCORE_COMPARISON_THRESHOLD",
    "GeneralMotifComparator",
    "MISSING_ADJUSTED_PVALUE_MESSAGE",
    "SUPPORTED_MOTIF_METRICS",
    "SUPPORTED_PROFILE_METRICS",
    "TomtomComparator",
    "UniversalMotifComparator",
    "compare",
    "compare_one_to_many",
    "comparison_column_for_criterion",
    "create_comparator_config",
    "default_threshold_for_criterion",
    "mimosa_compare_one_to_many",
]
