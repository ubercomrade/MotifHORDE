"""Compatibility facade for the generic MIMOSA motif runtime."""

from __future__ import annotations

import os
from typing import Optional

import numpy as np

from mimosa.handlers import pwm_model_from_pfm
from mimosa.models import (
    GenericModel,
    get_model_handler,
    read_model as _read_model,
    read_models,
    register_model_handler,
    registry,
    write_model,
)
from mimosa.scanning import (
    StrandMode,
    calculate_threshold_table as _calculate_threshold_table,
    get_frequencies,
    get_score_bounds,
    get_scores,
    scan_model,
    scan_model_strands,
)
from mimosa.sites import get_pfm as _get_pfm
from mimosa.sites import get_sites as _get_sites
from motifhorde.io import read_meme as _read_legacy_meme


def read_model(path: str | os.PathLike, model_type: str, **kwargs) -> GenericModel:
    """Read a model while preserving motifhorde's legacy PWM .txt support."""
    path_str = os.fspath(path)
    _, ext = os.path.splitext(path_str.lower())

    if model_type == "pwm" and ext in {".meme", ".txt"}:
        pfm, info, _ = _read_legacy_meme(path_str, index=kwargs.get("index", 0))
        name, length = info
        return pwm_model_from_pfm(pfm, name, int(length))

    return _read_model(path_str, model_type, **kwargs)


def calculate_threshold_table(
    model: GenericModel,
    sequences,
    strand: StrandMode = "best",
) -> np.ndarray:
    """Calculate motifhorde-compatible score-to-log-tail lookup values."""
    return _calculate_threshold_table(model, sequences, strand=strand)


def _resolve_legacy_site_strand(mode: str, strand: Optional[StrandMode]) -> StrandMode:
    """Keep old motifhorde defaults while allowing explicit MIMOSA strand modes."""
    if strand is not None:
        return strand
    if str(mode).lower() == "threshold":
        return "both"
    return "best"


def _resolve_legacy_threshold_table(
    model: GenericModel,
    sequences,
    mode: str,
    strand: Optional[StrandMode],
    background_sequences,
    threshold_table: Optional[np.ndarray],
) -> Optional[np.ndarray]:
    """Use the old best-strand calibration unless callers opt into a strand."""
    if threshold_table is not None or str(mode).lower() != "threshold":
        return threshold_table

    calibration_sequences = background_sequences if background_sequences is not None else sequences
    calibration_strand = "best" if strand is None else strand
    return calculate_threshold_table(model, calibration_sequences, strand=calibration_strand)


def get_sites(
    model: GenericModel,
    sequences,
    mode: str = "best",
    fpr_threshold: Optional[float] = None,
    background_sequences=None,
    threshold_table: Optional[np.ndarray] = None,
    *,
    strand: Optional[StrandMode] = None,
):
    """Find motif binding sites with motifhorde's historical defaults."""
    resolved_table = _resolve_legacy_threshold_table(
        model,
        sequences,
        mode,
        strand,
        background_sequences,
        threshold_table,
    )
    return _get_sites(
        model,
        sequences,
        mode=mode,
        fpr_threshold=fpr_threshold,
        strand=_resolve_legacy_site_strand(mode, strand),
        background_sequences=background_sequences,
        threshold_table=resolved_table,
    )


def get_pfm(
    model: GenericModel,
    sequences,
    mode: str = "best",
    fpr_threshold: Optional[float] = None,
    background_sequences=None,
    threshold_table: Optional[np.ndarray] = None,
    top_fraction: Optional[float] = None,
    pseudocount: float = 0.25,
    *,
    strand: Optional[StrandMode] = None,
) -> np.ndarray:
    """Construct a PFM from binding sites with motifhorde's historical defaults."""
    resolved_table = _resolve_legacy_threshold_table(
        model,
        sequences,
        mode,
        strand,
        background_sequences,
        threshold_table,
    )
    return _get_pfm(
        model,
        sequences,
        mode=mode,
        fpr_threshold=fpr_threshold,
        strand=_resolve_legacy_site_strand(mode, strand),
        background_sequences=background_sequences,
        threshold_table=resolved_table,
        top_fraction=top_fraction,
        pseudocount=pseudocount,
    )


__all__ = [
    "GenericModel",
    "StrandMode",
    "calculate_threshold_table",
    "get_frequencies",
    "get_model_handler",
    "get_pfm",
    "get_score_bounds",
    "get_scores",
    "get_sites",
    "read_model",
    "read_models",
    "register_model_handler",
    "registry",
    "scan_model",
    "scan_model_strands",
    "write_model",
]
