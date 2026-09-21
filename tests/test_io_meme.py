from __future__ import annotations

import re

import numpy as np
import pytest

from motifhorde.io import write_meme


def _pfm(nsites: float) -> np.ndarray:
    return np.full((4, 3), nsites / 4.0, dtype=np.float32)


def test_write_meme_reports_real_nsites(tmp_path):
    path = tmp_path / "export.meme"

    write_meme([_pfm(6.0)], [("Motif", 3, 6)], path)

    match = re.search(r"nsites= (\d+)", path.read_text())
    assert match is not None
    assert int(match.group(1)) == 6


def test_write_meme_rejects_divergent_column_sums(tmp_path):
    broken = _pfm(6.0)
    broken[:, 2] *= 2.0

    with pytest.raises(ValueError, match="column sums"):
        write_meme([broken], [("Motif", 3, 9)], tmp_path / "export.meme")
