from __future__ import annotations

import os
import subprocess

import numpy as np

from hordemotifs.discovery import MemeDiscoveryTool
from hordemotifs.io import write_meme
from hordemotifs.models import GenericModel


def _write_two_motif_meme(path):
    pfm = np.full((4, 4), 0.25, dtype=np.float32)
    write_meme([pfm, pfm], [("A", 4), ("B", 4)], path)


def test_meme_discovery_reads_generic_models(monkeypatch, tmp_path):
    calls = []

    def fake_run_checked(args, cwd=None):
        calls.append(args)
        meme_path = tmp_path / "source.meme"
        _write_two_motif_meme(meme_path)
        return subprocess.CompletedProcess(args, 0, stdout=meme_path.read_text(), stderr="")

    monkeypatch.setattr("hordemotifs.discovery.run_checked", fake_run_checked)

    tool = MemeDiscoveryTool(command="meme-bin", objfun="de", seed=7, threads=2)
    motifs = tool.discover("fg.fa", "bg.fa", os.fspath(tmp_path / "out"), 3, length=4)

    assert [motif.name for motif in motifs] == ["Meme-1", "Meme-2"]
    assert all(isinstance(motif, GenericModel) for motif in motifs)
    assert all(motif.type_key == "pwm" for motif in motifs)
    assert calls[0][:6] == ["meme-bin", "fg.fa", "-dna", "-revcomp", "-neg", "bg.fa"]
    assert "-seed" in calls[0]
    assert "-p" in calls[0]


def test_meme_discovery_missing_output_returns_empty(monkeypatch, tmp_path):
    def fake_run_checked(args, cwd=None):
        return subprocess.CompletedProcess(args, 0, stdout="", stderr="")

    monkeypatch.setattr("hordemotifs.discovery.run_checked", fake_run_checked)

    motifs = MemeDiscoveryTool(command="meme-bin").discover("fg.fa", "bg.fa", os.fspath(tmp_path), 2, length=4)

    assert motifs == []
