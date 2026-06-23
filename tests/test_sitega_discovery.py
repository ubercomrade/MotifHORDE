from __future__ import annotations

import os
import sys
from types import SimpleNamespace

import pytest

from motifhorde.discovery import SitegaDiscoveryTool


def _write_fasta(path) -> None:
    path.write_text(">seq\nACGT\n")


def _write_sitega_mat(path, length: int) -> None:
    path.write_text(
        "\n".join(
            [
                "SiteGATest",
                "1\tLPD count",
                f"{length}\tModel length",
                "0.0\tMinimum",
                "1.0\tRazmah",
                "0\t2\t1.5\t0\tac",
            ]
        )
    )


def test_sitega_discovery_uses_returned_mat_path(monkeypatch, tmp_path):
    foreground = tmp_path / "foreground.fa"
    background = tmp_path / "background.fa"
    output_dir = tmp_path / "sitega-output"
    returned_mat = tmp_path / "returned.mat"
    fallback_mat = output_dir / "train.fa_mat_fallback"
    calls = []

    _write_fasta(foreground)
    _write_fasta(background)

    def train(**kwargs):
        calls.append(kwargs)
        output_dir.mkdir(exist_ok=True)
        _write_sitega_mat(returned_mat, length=4)
        _write_sitega_mat(fallback_mat, length=6)
        return 0, os.fspath(returned_mat), "", 0.25

    monkeypatch.setitem(sys.modules, "sitega", SimpleNamespace(train=train))

    motifs = SitegaDiscoveryTool().discover(
        os.fspath(foreground),
        os.fspath(background),
        os.fspath(output_dir),
        number_of_motifs=5,
        length=4,
        lpd=10,
    )

    assert calls == [
        {
            "fg_path": os.fspath(foreground),
            "bg_path": os.fspath(background),
            "max_lpd": 6,
            "motif_len": 4,
            "size": 10,
            "olig_bg": 6,
            "infc": 1,
            "out_path": os.fspath(output_dir) + os.sep,
            "max_peak_len": 5000,
            "log_file": "sitega.log",
        }
    ]
    assert [(motif.name, motif.length) for motif in motifs] == [("Sitega-1", 4)]
    assert not (output_dir / "train.fa").exists()


def test_sitega_discovery_falls_back_to_glob_for_empty_mat_path(monkeypatch, tmp_path):
    foreground = tmp_path / "foreground.fa"
    background = tmp_path / "background.fa"
    output_dir = tmp_path / "sitega-output"
    fallback_mat = output_dir / "train.fa_mat0.mat"

    _write_fasta(foreground)
    _write_fasta(background)

    def train(**kwargs):
        output_dir.mkdir(exist_ok=True)
        _write_sitega_mat(fallback_mat, length=4)
        return 0, "", "", 0.25

    monkeypatch.setitem(sys.modules, "sitega", SimpleNamespace(train=train))

    motifs = SitegaDiscoveryTool().discover(
        os.fspath(foreground),
        os.fspath(background),
        os.fspath(output_dir),
        number_of_motifs=5,
        length=4,
        lpd=10,
    )

    assert [(motif.name, motif.length) for motif in motifs] == [("Sitega-1", 4)]


def test_sitega_discovery_rejects_failed_training(monkeypatch, tmp_path):
    foreground = tmp_path / "foreground.fa"
    background = tmp_path / "background.fa"
    output_dir = tmp_path / "sitega-output"

    _write_fasta(foreground)
    _write_fasta(background)

    def train(**kwargs):
        return 9, "", "", 0.0

    monkeypatch.setitem(sys.modules, "sitega", SimpleNamespace(train=train))

    with pytest.raises(RuntimeError, match="return code 9"):
        SitegaDiscoveryTool().discover(
            os.fspath(foreground),
            os.fspath(background),
            os.fspath(output_dir),
            number_of_motifs=5,
            length=4,
            lpd=10,
        )
