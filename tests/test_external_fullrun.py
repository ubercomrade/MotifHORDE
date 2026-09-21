from __future__ import annotations

import os
from pathlib import Path

import pytest

from motifhorde.discovery import (
    DimontDiscoveryTool,
    MemeDiscoveryTool,
    SlimDiscoveryTool,
)

pytestmark = [pytest.mark.external, pytest.mark.fullrun]


def _require_fullrun(command_env: str) -> str:
    if os.environ.get("HORDEMOTIFS_RUN_FULLRUN") != "1":
        pytest.skip("set HORDEMOTIFS_RUN_FULLRUN=1 to run external smoke tests")
    command = os.environ.get(command_env)
    if not command:
        pytest.skip(f"set {command_env} to run this external smoke test")
    return command


def _require_jar(jar_env: str, default: str) -> str:
    if os.environ.get("HORDEMOTIFS_RUN_FULLRUN") != "1":
        pytest.skip("set HORDEMOTIFS_RUN_FULLRUN=1 to run external smoke tests")
    jar = os.environ.get(jar_env, default)
    if not Path(jar).is_file():
        pytest.skip(f"Jstacs JAR is missing: {jar}")
    return jar


def _small_inputs() -> tuple[str, str]:
    root = Path(__file__).parent / "test_data" / "small_pipeline"
    return os.fspath(root / "foreground.fa"), os.fspath(root / "background.fa")


def test_meme_full_pipeline_smoke(tmp_path):
    command = _require_fullrun("HORDEMOTIFS_MEME_COMMAND")
    foreground, background = _small_inputs()
    models = MemeDiscoveryTool(command=command).discover(
        foreground,
        background,
        os.fspath(tmp_path),
        1,
        length=4,
    )
    assert len(models) == 1
    assert models[0].motif_length == 4


@pytest.mark.parametrize(
    ("tool_cls", "jar_env", "jar_default"),
    [
        pytest.param(
            DimontDiscoveryTool,
            "HORDEMOTIFS_DIMONT_JAR",
            "/home/anton/Programs/Jstacs/Dimont.jar",
            id="dimont-/home/anton/Programs/Jstacs/Dimont.jar",
        ),
        pytest.param(
            SlimDiscoveryTool,
            "HORDEMOTIFS_SLIM_JAR",
            "/home/anton/Programs/Jstacs/SlimDimont.jar",
            id="slim-/home/anton/Programs/Jstacs/SlimDimont.jar",
        ),
    ],
)
def test_jstacs_full_pipeline_smoke(tool_cls, jar_env, jar_default, tmp_path):
    jar = _require_jar(jar_env, jar_default)
    java = os.environ.get("HORDEMOTIFS_JAVA_COMMAND", "java")
    foreground, background = _small_inputs()
    models = tool_cls(jar_path=jar, java_command=java, threads=1).discover(
        foreground,
        background,
        os.fspath(tmp_path),
        1,
        length=4,
    )
    assert len(models) == 1
    assert models[0].motif_length == 4
