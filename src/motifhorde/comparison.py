"""Pipeline-facing comparison wrappers backed by MIMOSA."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

import pandas as pd
from mimosa.comparison import (
    SUPPORTED_MOTIF_METRICS,
    SUPPORTED_PROFILE_METRICS,
    compare,
    compare_one_to_many,
    create_comparator_config,
)
from mimosa.types import ComparatorConfig, ComparisonResult


class GeneralMotifComparator:
    """Small pipeline-facing comparator interface."""

    def __init__(self, name: str) -> None:
        self.name = name

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
        *,
        pvalue: bool = False,
        null_distribution: str | Path | dict[str, Any] | None = None,
        null_search_dirs: Iterable[str | Path] | None = None,
        effective_number_of_targets: int | None = None,
    ) -> None:
        super().__init__(name=f"TomtomComparator_{metric.upper()}")
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
                compare_one_to_many(
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
        filter_type: str | None = None,
        filter_threshold: float = 0.05,
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
        super().__init__(name)
        self.filter_type = filter_type
        self.filter_threshold = filter_threshold
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
                compare_one_to_many(
                    motif, motifs_2, "profile", self.config, sequences=sequences
                )
            )
        frame = _records_to_frame(records)
        if self.filter_type is None or frame.empty:
            return frame
        if self.filter_type not in frame.columns:
            raise ValueError(
                f"Comparison results do not contain filter column: {self.filter_type}"
            )
        if self.filter_type == "p-value":
            return frame[frame[self.filter_type] <= self.filter_threshold].reset_index(
                drop=True
            )
        return frame[frame[self.filter_type] >= self.filter_threshold].reset_index(
            drop=True
        )


def _records_to_frame(records):
    return pd.DataFrame.from_records(records)


__all__ = [
    "ComparatorConfig",
    "ComparisonResult",
    "GeneralMotifComparator",
    "SUPPORTED_MOTIF_METRICS",
    "SUPPORTED_PROFILE_METRICS",
    "TomtomComparator",
    "UniversalMotifComparator",
    "compare",
    "compare_one_to_many",
    "create_comparator_config",
]
