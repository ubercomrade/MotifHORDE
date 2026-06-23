from __future__ import annotations

import pytest

from motifhorde.cli import (
    check_dependencies,
    create_arg_parser,
    setup_comparator,
    setup_discovery_params,
    setup_discovery_tool,
    validate_args,
)
from motifhorde.comparison import TomtomComparator, UniversalMotifComparator
from motifhorde.discovery import (
    BammDiscoveryTool,
    DimontDiscoveryTool,
    MemeDiscoveryTool,
    SitegaDiscoveryTool,
    SlimDiscoveryTool,
    StremeDiscoveryTool,
)


def _parse(tool: str, *extra: str, validate: bool = True):
    parser = create_arg_parser()
    args = parser.parse_args(["fg.fa", "bg.fa", "prom.fa", "out", "-t", tool, *extra])
    if validate:
        validate_args(parser, args)
    return args


def test_cli_help_includes_new_tools_and_options():
    help_text = create_arg_parser().format_help()

    for text in [
        "meme",
        "dimont",
        "slim",
        "--meme-command",
        "--dimont-jar",
        "--slim-jar",
        "--java-xmx",
        "--mimosa-metric",
        "--comparison-criterion",
        "--mimosa-null-distribution",
        "--jobs",
    ]:
        assert text in help_text


def test_cli_help_excludes_removed_options():
    help_text = create_arg_parser().format_help()

    for text in [
        "continuous",
        "--c-metric",
        "--c-filter",
        "--c-threshold",
        "--c-search-range",
        "--c-jobs",
        "--tomtom-jobs",
        "--meme-p",
        "--jstacs-threads",
    ]:
        assert text not in help_text


def test_parser_accepts_mimosa_and_rejects_continuous():
    assert _parse("streme", "-c", "mimosa").comparator == "mimosa"

    parser = create_arg_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["fg.fa", "bg.fa", "prom.fa", "out", "-c", "continuous"])


def test_pvalue_comparison_requires_null_distribution():
    parser = create_arg_parser()
    args = parser.parse_args(
        [
            "fg.fa",
            "bg.fa",
            "prom.fa",
            "out",
            "-c",
            "mimosa",
            "--comparison-criterion",
            "p-value",
        ]
    )

    with pytest.raises(SystemExit):
        validate_args(parser, args)


@pytest.mark.parametrize("jobs", ["0", "-2"])
def test_jobs_validation_rejects_non_positive_values(jobs):
    parser = create_arg_parser()
    args = parser.parse_args(["fg.fa", "bg.fa", "prom.fa", "out", "--jobs", jobs])

    with pytest.raises(SystemExit):
        validate_args(parser, args)


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


@pytest.mark.parametrize(
    ("tool", "expected_attr"),
    [
        ("meme", "threads"),
        ("dimont", "threads"),
        ("slim", "threads"),
        ("sitega", "threads"),
    ],
)
def test_setup_discovery_tool_propagates_jobs(tool, expected_attr):
    discovery_tool = setup_discovery_tool(_parse(tool, "--jobs", "3"))

    assert getattr(discovery_tool, expected_attr) == 3


@pytest.mark.parametrize("tool", ["meme", "dimont", "slim", "sitega"])
def test_setup_discovery_tool_resolves_auto_jobs_for_external_tools(monkeypatch, tool):
    monkeypatch.setattr("motifhorde.cli.os.cpu_count", lambda: 7)

    discovery_tool = setup_discovery_tool(_parse(tool, "--jobs", "-1"))

    assert discovery_tool.threads == 7


def test_setup_comparator_configures_tomtom_jobs():
    comparator = setup_comparator(_parse("streme", "--jobs", "4"))

    assert isinstance(comparator, TomtomComparator)
    assert comparator.config.n_jobs == 4
    assert comparator.comparison_criterion == "score"
    assert comparator.comparison_threshold == 0.9


def test_setup_comparator_configures_mimosa_score_mode():
    comparator = setup_comparator(
        _parse(
            "streme",
            "-c",
            "mimosa",
            "--mimosa-metric",
            "dice",
            "--jobs",
            "4",
        )
    )

    assert isinstance(comparator, UniversalMotifComparator)
    assert comparator.config.metric == "dice"
    assert comparator.config.n_jobs == 4
    assert comparator.config.pvalue is False
    assert comparator.comparison_criterion == "score"
    assert comparator.comparison_threshold == 0.9


def test_setup_comparator_configures_mimosa_pvalue_mode():
    comparator = setup_comparator(
        _parse(
            "streme",
            "-c",
            "mimosa",
            "--comparison-criterion",
            "p-value",
            "--mimosa-null-distribution",
            "profile-null.joblib",
        )
    )

    assert isinstance(comparator, UniversalMotifComparator)
    assert comparator.config.pvalue is True
    assert comparator.config.null_distribution == "profile-null.joblib"
    assert comparator.comparison_criterion == "p-value"
    assert comparator.comparison_column == "adj.p-value"
    assert comparator.comparison_threshold == 0.05


def test_setup_discovery_params_keeps_tool_specific_values():
    bamm = setup_discovery_params(_parse("bamm", "-l", "8", "-o", "1,2"))
    sitega = setup_discovery_params(_parse("sitega", "-l", "8", "--lpd", "10"))
    meme = setup_discovery_params(_parse("meme", "-l", "8"))

    assert bamm == {"length": [8], "order": [1, 2]}
    assert sitega == {"length": [8], "lpd": [10]}
    assert meme == {"length": [8]}


def test_sitega_dependency_check_does_not_require_executable(monkeypatch):
    def fail_which(command: str) -> str | None:
        raise AssertionError(f"unexpected PATH lookup for {command}")

    monkeypatch.setattr("motifhorde.cli.shutil.which", fail_which)

    check_dependencies(_parse("sitega"))
