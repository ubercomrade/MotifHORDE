"""Small helpers for resolving and running external commands."""

from __future__ import annotations

import os
import shutil
import subprocess

DEFAULT_MEME_COMMAND = "meme"
DEFAULT_STREME_COMMAND = "streme"
PACKAGE_BIN_DIR = os.path.join(os.path.dirname(__file__), "bin")
DEFAULT_DIMONT_JAR = os.path.join(PACKAGE_BIN_DIR, "Dimont.jar")
DEFAULT_SLIM_JAR = os.path.join(PACKAGE_BIN_DIR, "SlimDimont.jar")
DEFAULT_BAMM_COMMAND = "BaMMmotif"


def resolve_command(
    command_name: str,
    fallback_path: str | None = None,
    env_var: str | None = None,
) -> str:
    """Resolve a command from environment, PATH, then fallback."""
    if env_var:
        env_value = os.environ.get(env_var)
        if env_value:
            return env_value

    path_value = shutil.which(command_name)
    if path_value:
        return path_value

    return fallback_path or command_name


def resolve_existing_path(
    cli_value: str | None,
    env_var: str,
    default_path: str,
    label: str,
) -> str:
    """Resolve a file path and require that it exists."""
    path = cli_value or os.environ.get(env_var) or default_path
    if not os.path.exists(path):
        raise FileNotFoundError(f"{label} not found: {path}")
    return path


def run_checked(args: list[str], cwd: str | None = None) -> subprocess.CompletedProcess:
    """Run a subprocess and raise a detailed error on failure."""
    result = subprocess.run(args, cwd=cwd, shell=False, capture_output=True, text=True)
    if result.returncode != 0:
        command = " ".join(args)
        raise RuntimeError(
            f"Command failed: {command}\n"
            f"Return code: {result.returncode}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )
    return result
