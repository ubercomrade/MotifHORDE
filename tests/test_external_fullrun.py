from __future__ import annotations

import os
import shutil
import subprocess

import pytest

from hordemotifs.external import (
    DEFAULT_DIMONT_JAR,
    DEFAULT_MEME_COMMAND,
    DEFAULT_SLIM_JAR,
    resolve_command,
)

SMALL_DATA = "tests/test_data/small_pipeline"
JSTACS_EXAMPLE = "/home/anton/Programs/Jstacs/dimont-example.fa"


def _run_cli(tool: str, output_dir, foreground: str, background: str, promoters: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            "uv",
            "run",
            "hordeMotifs",
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
            "--tomtom-perm",
            "0",
            "--jstacs-threads",
            "1",
            "--dimont-starts",
            "1",
            "--slim-starts",
            "1",
        ],
        shell=False,
        capture_output=True,
        text=True,
    )


@pytest.mark.external
@pytest.mark.meme
@pytest.mark.fullrun
def test_meme_full_pipeline_smoke(tmp_path):
    if os.environ.get("HORDEMOTIFS_RUN_FULLRUN") != "1":
        pytest.skip("Set HORDEMOTIFS_RUN_FULLRUN=1 to run external full pipeline smoke tests")
    command = resolve_command("meme", DEFAULT_MEME_COMMAND, "HORDEMOTIFS_MEME_COMMAND")
    if not (os.path.exists(command) or shutil.which(command)):
        pytest.skip("MEME executable is not available")

    result = _run_cli(
        "meme",
        tmp_path / "meme-out",
        f"{SMALL_DATA}/foreground.fa",
        f"{SMALL_DATA}/background.fa",
        f"{SMALL_DATA}/promoters.fa",
    )

    assert result.returncode == 0, result.stderr


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
        pytest.skip("Set HORDEMOTIFS_RUN_FULLRUN=1 to run external full pipeline smoke tests")
    if shutil.which("java") is None or not os.path.exists(jar_path):
        pytest.skip(f"{tool} dependencies are not available")
    if not os.path.exists(JSTACS_EXAMPLE):
        pytest.skip("Jstacs example FASTA is not available")

    result = _run_cli(tool, tmp_path / f"{tool}-out", JSTACS_EXAMPLE, JSTACS_EXAMPLE, JSTACS_EXAMPLE)

    assert result.returncode == 0, result.stderr
