from __future__ import annotations

import os
import subprocess

from motifhorde.discovery import DimontDiscoveryTool, _build_dimont_args
from motifhorde.models import GenericModel


def test_build_dimont_args_contains_jstacs_parameters():
    args = _build_dimont_args(
        "java",
        "2G",
        "Dimont.jar",
        "/tmp/out",
        "train.annot.fa",
        8,
        "position",
        "value",
        -1,
        0,
        4.0,
        20,
        3,
    )

    assert args[:5] == ["java", "-Djava.awt.headless=true", "-Xmx2G", "-jar", "Dimont.jar"]
    assert "home=/tmp/out" in args
    assert "data=train.annot.fa" in args
    assert "infix=dimont" in args
    assert "motifWidth=8" in args
    assert "threads=3" in args


def test_dimont_discovery_writes_fasta_and_reads_generic_model(monkeypatch, tmp_path):
    jar_path = tmp_path / "Dimont.jar"
    jar_path.write_text("")
    foreground = tmp_path / "fg.fa"
    foreground.write_text(">x\nACGTACGT\n")
    output_dir = tmp_path / "out"

    def fake_run_checked(args, cwd=None):
        output_dir.mkdir(exist_ok=True)
        (output_dir / "result_dimont.xml").write_text("<xml />")
        return subprocess.CompletedProcess(args, 0, stdout="", stderr="")

    def fake_read_model(path, model_type, **kwargs):
        return GenericModel(model_type, "raw", [[[0.0]]], 1, {"kmer": 1})

    monkeypatch.setattr("motifhorde.discovery.run_checked", fake_run_checked)
    monkeypatch.setattr("motifhorde.discovery.read_model", fake_read_model)

    motifs = DimontDiscoveryTool(jar_path=os.fspath(jar_path)).discover(
        os.fspath(foreground),
        "bg.fa",
        os.fspath(output_dir),
        1,
        length=8,
    )

    assert (output_dir / "train.annot.fa").read_text() == "> position: 4; value: 1.0\nACGTACGT\n"
    assert len(motifs) == 1
    assert motifs[0].type_key == "dimont"
    assert motifs[0].name == "Dimont-1"
