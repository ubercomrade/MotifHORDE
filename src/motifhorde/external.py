"""Small helpers for resolving and running external commands."""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess

DEFAULT_MEME_COMMAND = "meme"
DEFAULT_STREME_COMMAND = "streme"
DEFAULT_BAMM_COMMAND = "BaMMmotif"
DEFAULT_SITEGA_COMMAND = "sitega"


def resolve_command(
    command_name: str,
    explicit: str | None = None,
    env_var: str | None = None,
) -> str:
    """Resolve a command as CLI value, environment value, then PATH."""
    if explicit:
        return explicit
    if env_var:
        env_value = os.environ.get(env_var)
        if env_value:
            return env_value

    path_value = shutil.which(command_name)
    if path_value:
        return path_value

    return command_name


def resolve_existing_path(
    cli_value: str | None,
    env_var: str,
    label: str,
) -> str:
    """Resolve a file path and require that it exists."""
    path = cli_value or os.environ.get(env_var)
    if not path:
        raise FileNotFoundError(f"{label} requires a CLI path or {env_var}")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"{label} not found: {path}")
    return path


def run_checked(args: list[str], cwd: str | None = None) -> subprocess.CompletedProcess:
    """Run a subprocess and raise a detailed error on failure."""
    try:
        result = subprocess.run(
            args, cwd=cwd, shell=False, capture_output=True, text=True
        )
    except OSError as exc:
        raise RuntimeError(
            f"Command failed to start: {shlex.join(args)}\nstdout:\n\nstderr:\n{exc}"
        ) from exc
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed: {shlex.join(args)}\n"
            f"Return code: {result.returncode}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )
    return result
