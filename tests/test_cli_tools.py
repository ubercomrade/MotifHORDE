from __future__ import annotations

import pytest

from motifhorde.cli import create_arg_parser, setup_discovery_params, setup_discovery_tool
from motifhorde.discovery import (
    BammDiscoveryTool,
    DimontDiscoveryTool,
    MemeDiscoveryTool,
    SitegaDiscoveryTool,
    SlimDiscoveryTool,
    StremeDiscoveryTool,
)


def _parse(tool: str, *extra: str):
    parser = create_arg_parser()
    return parser.parse_args(["fg.fa", "bg.fa", "prom.fa", "out", "-t", tool, *extra])


def test_cli_help_includes_new_tools_and_options():
    help_text = create_arg_parser().format_help()

    for text in ["meme", "dimont", "slim", "--meme-command", "--dimont-jar", "--slim-jar", "--java-xmx"]:
        assert text in help_text


@pytest.mark.parametrize(
    ("tool", "expected_type"),
    [
        ("streme", StremeDiscoveryTool),
        ("meme", MemeDiscoveryTool),
        ("bamm", BammDiscoveryTool),
        ("dimont", DimontDiscoveryTool),
        ("slim", SlimDiscoveryTool),
        ("sitega", SitegaDiscoveryTool),
    ],
)
def test_setup_discovery_tool_returns_selected_tool(tool, expected_type):
    assert isinstance(setup_discovery_tool(_parse(tool)), expected_type)


def test_setup_discovery_params_keeps_tool_specific_values():
    bamm = setup_discovery_params(_parse("bamm", "-l", "8", "-o", "1,2"))
    sitega = setup_discovery_params(_parse("sitega", "-l", "8", "--lpd", "10"))
    meme = setup_discovery_params(_parse("meme", "-l", "8"))

    assert bamm == {"length": [8], "order": [1, 2]}
    assert sitega == {"length": [8], "lpd": [10]}
    assert meme == {"length": [8]}
