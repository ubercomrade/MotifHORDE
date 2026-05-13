from __future__ import annotations

import os
import subprocess

import numpy as np

from motifhorde.discovery import MemeDiscoveryTool, StremeDiscoveryTool
from motifhorde.io import write_meme
from motifhorde.models import GenericModel


def _write_motifs_meme(path, lengths):
    pfms = [np.full((4, length), 0.25, dtype=np.float32) for length in lengths]
    info = [(f"motif-{index}", length) for index, length in enumerate(lengths, start=1)]
    write_meme(pfms, info, path)


def test_meme_discovery_reads_generic_models(monkeypatch, tmp_path):
    calls = []

    def fake_run_checked(args, cwd=None):
        calls.append(args)
        meme_path = tmp_path / "source.meme"
        _write_motifs_meme(meme_path, [4, 4])
        return subprocess.CompletedProcess(args, 0, stdout=meme_path.read_text(), stderr="")

    monkeypatch.setattr("motifhorde.discovery.run_checked", fake_run_checked)

    tool = MemeDiscoveryTool(command="meme-bin", objfun="de", seed=7, threads=2)
    motifs = tool.discover("fg.fa", "bg.fa", os.fspath(tmp_path / "out"), 3, length=4)

    assert [motif.name for motif in motifs] == ["Meme-1", "Meme-2"]
    assert all(isinstance(motif, GenericModel) for motif in motifs)
    assert all(motif.type_key == "pwm" for motif in motifs)
    assert calls[0][:6] == ["meme-bin", "fg.fa", "-dna", "-revcomp", "-neg", "bg.fa"]
    assert "-minw" in calls[0]
    assert "-maxw" in calls[0]
    assert "-nomatrim" in calls[0]
    assert "-seed" in calls[0]
    assert "-p" in calls[0]


def test_meme_discovery_missing_output_returns_empty(monkeypatch, tmp_path):
    def fake_run_checked(args, cwd=None):
        return subprocess.CompletedProcess(args, 0, stdout="", stderr="")

    monkeypatch.setattr("motifhorde.discovery.run_checked", fake_run_checked)

    motifs = MemeDiscoveryTool(command="meme-bin").discover("fg.fa", "bg.fa", os.fspath(tmp_path), 2, length=4)

    assert motifs == []


def test_streme_discovery_sets_strict_width_bounds(monkeypatch, tmp_path):
    calls = []

    def fake_run_checked(args, cwd=None):
        calls.append(args)
        meme_path = tmp_path / "source.meme"
        _write_motifs_meme(meme_path, [4])
        return subprocess.CompletedProcess(args, 0, stdout=meme_path.read_text(), stderr="")

    monkeypatch.setattr("motifhorde.discovery.run_checked", fake_run_checked)

    motifs = StremeDiscoveryTool(command="streme-bin").discover(
        "fg.fa",
        "bg.fa",
        os.fspath(tmp_path / "out"),
        1,
        length=4,
    )

    assert [motif.name for motif in motifs] == ["Streme-1"]
    assert "--minw" in calls[0]
    assert "--maxw" in calls[0]
    assert "--w" not in calls[0]


def test_discovery_filters_motifs_with_unexpected_length(monkeypatch, tmp_path):
    def fake_run_checked(args, cwd=None):
        meme_path = tmp_path / "source.meme"
        _write_motifs_meme(meme_path, [4, 5])
        return subprocess.CompletedProcess(args, 0, stdout=meme_path.read_text(), stderr="")

    monkeypatch.setattr("motifhorde.discovery.run_checked", fake_run_checked)

    motifs = StremeDiscoveryTool(command="streme-bin").discover(
        "fg.fa",
        "bg.fa",
        os.fspath(tmp_path / "out"),
        2,
        length=4,
    )

    assert [motif.name for motif in motifs] == ["Streme-1"]
    assert [motif.length for motif in motifs] == [4]
