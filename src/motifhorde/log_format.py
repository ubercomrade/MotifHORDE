"""Small helpers for stable structured log fields."""

from __future__ import annotations

from typing import Any


def format_log_params(params: dict[str, Any]) -> str:
    return ",".join(f"{key}:{params[key]}" for key in sorted(params))


def format_elapsed(seconds: float) -> str:
    return f"{seconds:.1f}s"
