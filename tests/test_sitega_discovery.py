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
            "num_threads": 0,
            "seed": 0,
            "pop_size": 100,
            "num_motifs": 20,
        }
    ]
    assert [(motif.name, motif.length) for motif in motifs] == [("Sitega-1", 4)]
    assert not (output_dir / "train.fa").exists()


def test_sitega_discovery_normalizes_returned_path_without_extension(
    monkeypatch, tmp_path
):
    foreground = tmp_path / "foreground.fa"
    background = tmp_path / "background.fa"
    output_dir = tmp_path / "sitega-output"
    returned_mat = tmp_path / "returned"

    _write_fasta(foreground)
    _write_fasta(background)

    def train(**kwargs):
        output_dir.mkdir(exist_ok=True)
        _write_sitega_mat(returned_mat, length=4)
        return 0, os.fspath(returned_mat), "", 0.25

    monkeypatch.setitem(sys.modules, "sitega", SimpleNamespace(train=train))

    motifs = SitegaDiscoveryTool().discover(
        os.fspath(foreground),
        os.fspath(background),
        os.fspath(output_dir),
        number_of_motifs=1,
        length=4,
        lpd=10,
    )

    assert [(motif.name, motif.length) for motif in motifs] == [("Sitega-1", 4)]
    assert not returned_mat.exists()
    assert returned_mat.with_suffix(".mat").exists()


def test_sitega_discovery_passes_thread_count(monkeypatch, tmp_path):
    foreground = tmp_path / "foreground.fa"
    background = tmp_path / "background.fa"
    output_dir = tmp_path / "sitega-output"
    returned_mat = tmp_path / "returned.mat"
    calls = []

    _write_fasta(foreground)
    _write_fasta(background)

    def train(**kwargs):
        calls.append(kwargs)
        output_dir.mkdir(exist_ok=True)
        _write_sitega_mat(returned_mat, length=4)
        return 0, os.fspath(returned_mat), "", 0.25

    monkeypatch.setitem(sys.modules, "sitega", SimpleNamespace(train=train))

    SitegaDiscoveryTool(threads=3).discover(
        os.fspath(foreground),
        os.fspath(background),
        os.fspath(output_dir),
        number_of_motifs=1,
        length=4,
        lpd=10,
    )

    assert calls[0]["num_threads"] == 3


def test_sitega_discovery_passes_constructor_seed(monkeypatch, tmp_path):
    foreground = tmp_path / "foreground.fa"
    background = tmp_path / "background.fa"
    output_dir = tmp_path / "sitega-output"
    returned_mat = tmp_path / "returned.mat"
    calls = []

    _write_fasta(foreground)
    _write_fasta(background)

    def train(**kwargs):
        calls.append(kwargs)
        output_dir.mkdir(exist_ok=True)
        _write_sitega_mat(returned_mat, length=4)
        return 0, os.fspath(returned_mat), "", 0.25

    monkeypatch.setitem(sys.modules, "sitega", SimpleNamespace(train=train))

    SitegaDiscoveryTool(seed=123).discover(
        os.fspath(foreground),
        os.fspath(background),
        os.fspath(output_dir),
        number_of_motifs=1,
        length=4,
        lpd=10,
    )

    assert calls[0]["seed"] == 123


def test_sitega_discovery_kwargs_seed_overrides_constructor_seed(monkeypatch, tmp_path):
    foreground = tmp_path / "foreground.fa"
    background = tmp_path / "background.fa"
    output_dir = tmp_path / "sitega-output"
    returned_mat = tmp_path / "returned.mat"
    calls = []

    _write_fasta(foreground)
    _write_fasta(background)

    def train(**kwargs):
        calls.append(kwargs)
        output_dir.mkdir(exist_ok=True)
        _write_sitega_mat(returned_mat, length=4)
        return 0, os.fspath(returned_mat), "", 0.25

    monkeypatch.setitem(sys.modules, "sitega", SimpleNamespace(train=train))

    SitegaDiscoveryTool(seed=123).discover(
        os.fspath(foreground),
        os.fspath(background),
        os.fspath(output_dir),
        number_of_motifs=1,
        length=4,
        lpd=10,
        seed=456,
    )

    assert calls[0]["seed"] == 456


def test_sitega_discovery_real_backend_smoke(tmp_path):
    pytest.importorskip("sitega")

    output_dir = tmp_path / "sitega-output"
    foreground = "tests/test_data/small_pipeline/foreground.fa"
    background = "tests/test_data/small_pipeline/background.fa"

    motifs = SitegaDiscoveryTool(seed=123).discover(
        foreground,
        background,
        os.fspath(output_dir),
        number_of_motifs=1,
        length=8,
        lpd=10,
    )

    assert [(motif.name, motif.length) for motif in motifs] == [("Sitega-1", 8)]
    assert (output_dir / "foreground_mat1.mat").exists()
    assert "Pipeline completed successfully!" in (output_dir / "sitega.log").read_text()


def test_sitega_discovery_falls_back_to_glob_for_empty_mat_path(monkeypatch, tmp_path):
    foreground = tmp_path / "foreground.fa"
    background = tmp_path / "background.fa"
    output_dir = tmp_path / "sitega-output"
    fallback_mat = output_dir / "train_mat1"

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
    assert not fallback_mat.exists()
    assert fallback_mat.with_suffix(".mat").exists()


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
