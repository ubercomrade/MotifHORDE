from __future__ import annotations

import os
import subprocess

from motifhorde.discovery import SlimDiscoveryTool, _build_slim_args
from motifhorde.models import GenericModel


def test_build_slim_args_contains_jstacs_parameters():
    args = _build_slim_args(
        "java",
        "2G",
        "SlimDimont.jar",
        "/tmp/out",
        "train.annot.fa",
        8,
        "position",
        "value",
        -1,
        -5,
        False,
        2,
        None,
    )

    assert args[:5] == ["java", "-Djava.awt.headless=true", "-Xmx2G", "-jar", "SlimDimont.jar"]
    assert "infix=slim" in args
    assert "motifWidth=8" in args
    assert "motifOrder=-5" in args
    assert "starts=2" in args
    assert "modify=false" in args
    assert not any(arg.startswith("threads=") for arg in args)


def test_slim_discovery_writes_fasta_and_reads_generic_model(monkeypatch, tmp_path):
    jar_path = tmp_path / "SlimDimont.jar"
    jar_path.write_text("")
    foreground = tmp_path / "fg.fa"
    foreground.write_text(">x\nACGTACGT\n")
    output_dir = tmp_path / "out"

    def fake_run_checked(args, cwd=None):
        output_dir.mkdir(exist_ok=True)
        (output_dir / "result_slim.xml").write_text("<xml />")
        return subprocess.CompletedProcess(args, 0, stdout="", stderr="")

    def fake_read_model(path, model_type, **kwargs):
        return GenericModel(model_type, "raw", [[[0.0] * 8]], 8, {"kmer": 1})

    monkeypatch.setattr("motifhorde.discovery.run_checked", fake_run_checked)
    monkeypatch.setattr("motifhorde.discovery.read_model", fake_read_model)

    motifs = SlimDiscoveryTool(jar_path=os.fspath(jar_path)).discover(
        os.fspath(foreground),
        "bg.fa",
        os.fspath(output_dir),
        1,
        length=8,
    )

    assert (output_dir / "train.annot.fa").read_text() == "> position: 4; value: 1.0\nACGTACGT\n"
    assert len(motifs) == 1
    assert motifs[0].type_key == "slim"
    assert motifs[0].name == "Slim-1"
