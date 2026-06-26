from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Sequence

import pytest

from motifhorde.external import (
    DEFAULT_DIMONT_JAR,
    DEFAULT_MEME_COMMAND,
    DEFAULT_SLIM_JAR,
    DEFAULT_STREME_COMMAND,
    resolve_command,
)
from motifhorde.models import read_model

SMALL_DATA = "tests/test_data/small_pipeline"
REAL_FOREGROUND = "tests/test_data/PEAKS035260_GATA2_O09100_MACS2.fa"
REAL_BACKGROUND = "tests/test_data/PEAKS035260_GATA2_O09100_MACS2_gb.fa"
JSTACS_EXAMPLE = "/home/anton/Programs/Jstacs/dimont-example.fa"
REAL_DATA_RECORDS = 400


def _run_cli(
    tool: str,
    output_dir,
    foreground: str,
    background: str,
    promoters: str,
    extra_args: Sequence[str] = (),
) -> subprocess.CompletedProcess:
    python_command = os.environ.get("HORDEMOTIFS_TEST_PYTHON", sys.executable)
    return subprocess.run(
        [
            python_command,
            "-m",
            "motifhorde.cli",
            foreground,
            background,
            promoters,
            os.fspath(output_dir),
            "-t",
            tool,
            "-l",
            "6",
            "-n",
            "1",
            "--jobs",
            "1",
            "--dimont-starts",
            "1",
            "--slim-starts",
            "1",
            *extra_args,
        ],
        shell=False,
        capture_output=True,
        text=True,
    )


def _write_fasta_subset(source: str, target: Path, record_limit: int) -> None:
    records = 0
    with open(source, "r") as input_handle, open(target, "w") as output_handle:
        for line in input_handle:
            if line.startswith(">"):
                if records >= record_limit:
                    break
                records += 1
            if records:
                output_handle.write(line)

    if records < record_limit:
        raise ValueError(
            f"{source} contains {records} records, expected {record_limit}"
        )


def _write_run_log(log_path: Path, result: subprocess.CompletedProcess) -> None:
    log_path.write_text(
        "\n".join(
            [
                f"command: {' '.join(map(os.fspath, result.args))}",
                f"returncode: {result.returncode}",
                "",
                "stdout:",
                result.stdout,
                "",
                "stderr:",
                result.stderr,
            ]
        )
    )


def _prepare_real_data_subset(tmp_path: Path) -> tuple[str, str, str]:
    input_dir = tmp_path / "real-data"
    input_dir.mkdir()
    foreground = input_dir / "foreground.fa"
    background = input_dir / "background.fa"
    promoters = input_dir / "promoters.fa"

    _write_fasta_subset(REAL_FOREGROUND, foreground, REAL_DATA_RECORDS)
    _write_fasta_subset(REAL_BACKGROUND, background, REAL_DATA_RECORDS)
    _write_fasta_subset(REAL_BACKGROUND, promoters, REAL_DATA_RECORDS)

    return os.fspath(foreground), os.fspath(background), os.fspath(promoters)


def _java_command() -> str | None:
    env_value = os.environ.get("HORDEMOTIFS_JAVA_COMMAND")
    if env_value:
        return env_value

    path_value = shutil.which("java")
    if path_value:
        return path_value

    conda_exe = os.environ.get("CONDA_EXE")
    if conda_exe:
        conda_java = Path(conda_exe).with_name("java")
        if conda_java.exists():
            return os.fspath(conda_java)

    executable = Path(sys.executable).resolve()
    if "envs" in executable.parts:
        envs_index = executable.parts.index("envs")
        conda_root = Path(*executable.parts[:envs_index])
        conda_java = conda_root / "bin" / "java"
        if conda_java.exists():
            return os.fspath(conda_java)

    return None


def _skip_missing_real_data_dependency(tool: str) -> None:
    if tool == "streme":
        command = resolve_command(
            "streme", DEFAULT_STREME_COMMAND, "HORDEMOTIFS_STREME_COMMAND"
        )
        if not (os.path.exists(command) or shutil.which(command)):
            pytest.skip(f"STREME executable is not available: {command}")
    elif tool == "slim":
        if _java_command() is None or not os.path.exists(DEFAULT_SLIM_JAR):
            pytest.skip("SlimDimont dependencies are not available")
    elif tool == "sitega":
        pytest.importorskip("sitega")
    else:
        raise ValueError(f"Unsupported real-data tool: {tool}")


@pytest.mark.external
@pytest.mark.meme
@pytest.mark.fullrun
def test_meme_full_pipeline_smoke(tmp_path):
    if os.environ.get("HORDEMOTIFS_RUN_FULLRUN") != "1":
        pytest.skip(
            "Set HORDEMOTIFS_RUN_FULLRUN=1 to run external full pipeline smoke tests"
        )
    command = resolve_command("meme", DEFAULT_MEME_COMMAND, "HORDEMOTIFS_MEME_COMMAND")
    if not (os.path.exists(command) or shutil.which(command)):
        pytest.skip("MEME executable is not available")

    output_dir = tmp_path / "meme-out"
    result = _run_cli(
        "meme",
        output_dir,
        f"{SMALL_DATA}/foreground.fa",
        f"{SMALL_DATA}/background.fa",
        f"{SMALL_DATA}/promoters.fa",
    )

    assert result.returncode == 0, result.stderr
    models_dir = output_dir / "meme" / "motifs" / "models"
    model_paths = sorted(models_dir.glob("*.pkl"))
    assert model_paths
    assert (models_dir / "all_motifs_in_pfm_form.meme").exists()

    pkl_model = read_model(os.fspath(model_paths[0]), "pwm")
    meme_model = read_model(
        os.fspath(models_dir / "all_motifs_in_pfm_form.meme"), "pwm"
    )
    assert pkl_model.type_key == "pwm"
    assert meme_model.type_key == "pwm"


@pytest.mark.external
@pytest.mark.fullrun
@pytest.mark.parametrize("tool", ["streme", "slim", "sitega"])
def test_real_data_full_pipeline_smoke_with_verbose_log(tool, tmp_path):
    if os.environ.get("HORDEMOTIFS_RUN_FULLRUN") != "1":
        pytest.skip(
            "Set HORDEMOTIFS_RUN_FULLRUN=1 to run external full pipeline smoke tests"
        )
    _skip_missing_real_data_dependency(tool)

    foreground, background, promoters = _prepare_real_data_subset(tmp_path)
    output_dir = tmp_path / f"{tool}-real-out"
    log_dir = Path(
        os.environ.get("HORDEMOTIFS_FULLRUN_LOG_DIR", os.fspath(tmp_path / "logs"))
    )
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{tool}-real-data-fullrun.log"

    extra_args = [
        "--verbose",
        "--comparison-threshold",
        "0.0",
        "--lpd",
        "10",
    ]
    if tool == "slim":
        extra_args.extend(["--java-command", _java_command() or "java"])

    result = _run_cli(tool, output_dir, foreground, background, promoters, extra_args)
    _write_run_log(log_path, result)

    log_text = log_path.read_text()
    assert result.returncode == 0, log_text
    assert log_path.exists()
    assert "[pipeline] start" in log_text

    models_dir = output_dir / tool / "motifs" / "models"
    model_paths = sorted(models_dir.glob("*.pkl"))
    assert model_paths
    assert (models_dir / "all_motifs_in_pfm_form.meme").exists()


@pytest.mark.external
@pytest.mark.jstacs
@pytest.mark.fullrun
@pytest.mark.parametrize(
    ("tool", "jar_path"),
    [
        ("dimont", DEFAULT_DIMONT_JAR),
        ("slim", DEFAULT_SLIM_JAR),
    ],
)
def test_jstacs_full_pipeline_smoke(tool, jar_path, tmp_path):
    if os.environ.get("HORDEMOTIFS_RUN_FULLRUN") != "1":
        pytest.skip(
            "Set HORDEMOTIFS_RUN_FULLRUN=1 to run external full pipeline smoke tests"
        )
    if shutil.which("java") is None or not os.path.exists(jar_path):
        pytest.skip(f"{tool} dependencies are not available")
    if not os.path.exists(JSTACS_EXAMPLE):
        pytest.skip("Jstacs example FASTA is not available")

    result = _run_cli(
        tool, tmp_path / f"{tool}-out", JSTACS_EXAMPLE, JSTACS_EXAMPLE, JSTACS_EXAMPLE
    )

    assert result.returncode == 0, result.stderr
