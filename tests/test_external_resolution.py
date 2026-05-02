from __future__ import annotations

import os
import subprocess

import pytest

from motifhorde.external import resolve_command, resolve_existing_path, run_checked


def test_resolve_command_prefers_env(monkeypatch, tmp_path):
    env_command = tmp_path / "env-tool"
    monkeypatch.setenv("HORDEMOTIFS_TEST_COMMAND", os.fspath(env_command))

    assert resolve_command("missing-tool", "/fallback/tool", "HORDEMOTIFS_TEST_COMMAND") == os.fspath(env_command)


def test_resolve_command_uses_path_before_fallback(monkeypatch, tmp_path):
    tool = tmp_path / "demo-tool"
    tool.write_text("#!/bin/sh\n")
    tool.chmod(0o755)
    monkeypatch.setenv("PATH", os.fspath(tmp_path))

    assert resolve_command("demo-tool", "/fallback/tool") == os.fspath(tool)


def test_resolve_existing_path_order(monkeypatch, tmp_path):
    default_path = tmp_path / "default.jar"
    env_path = tmp_path / "env.jar"
    cli_path = tmp_path / "cli.jar"
    for path in [default_path, env_path, cli_path]:
        path.write_text("")

    monkeypatch.setenv("HORDEMOTIFS_TEST_JAR", os.fspath(env_path))

    assert resolve_existing_path(None, "HORDEMOTIFS_TEST_JAR", os.fspath(default_path), "test jar") == os.fspath(env_path)
    assert resolve_existing_path(os.fspath(cli_path), "HORDEMOTIFS_TEST_JAR", os.fspath(default_path), "test jar") == os.fspath(cli_path)


def test_resolve_existing_path_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="test jar not found"):
        resolve_existing_path(None, "HORDEMOTIFS_MISSING_JAR", os.fspath(tmp_path / "missing.jar"), "test jar")


def test_run_checked_error_contains_process_output():
    with pytest.raises(RuntimeError, match="Return code: 3") as exc_info:
        run_checked(["python", "-c", "import sys; print('out'); print('err', file=sys.stderr); sys.exit(3)"])

    message = str(exc_info.value)
    assert "out" in message
    assert "err" in message


def test_run_checked_success():
    result = run_checked(["python", "-c", "print('ok')"])
    assert isinstance(result, subprocess.CompletedProcess)
    assert result.stdout.strip() == "ok"
