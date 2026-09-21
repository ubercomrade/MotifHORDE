from __future__ import annotations

import pytest

from motifhorde.cli import (
    check_dependencies,
    create_arg_parser,
    main_cli,
    parse_range,
    setup_discovery_params,
    validate_args,
)


def _parse(argv: list[str]):
    return create_arg_parser().parse_args(["fg.fa", "bg.fa", "pr.fa", "out", *argv])


def test_cli_keeps_four_positionals_and_new_sitega_parameters():
    parser = create_arg_parser()
    args = parser.parse_args(
        [
            "foreground.fa",
            "background.fa",
            "promoters.fa",
            "out",
            "-t",
            "sitega",
            "--lpd",
            "4",
            "--features",
            "10,20",
        ]
    )

    assert setup_discovery_params(args) == {
        "length": [8, 12, 16, 20],
        "lpd": [4],
        "features": [10, 20],
    }
    help_text = parser.format_help()
    assert "--tomtom-metric" not in help_text
    assert "--pfm-mode" not in help_text
    assert "--sitega-command" in help_text


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("8", [8]),
        ("8,12,16", [8, 12, 16]),
        ("8-20-4", [8, 12, 16, 20]),
    ],
)
def test_parse_range_accepts_single_list_and_stepped_range(value, expected):
    assert parse_range(value) == expected


@pytest.mark.parametrize("value", ["8-20-0", "10-5-2", "20-8--4"])
def test_parse_range_rejects_zero_step_and_empty_results(value):
    with pytest.raises(ValueError):
        parse_range(value)


def _expect_system_exit(argv: list[str]) -> None:
    parser = create_arg_parser()
    args = parser.parse_args(["fg.fa", "bg.fa", "pr.fa", "out", *argv])
    with pytest.raises(SystemExit):
        validate_args(parser, args)


@pytest.mark.parametrize("jobs", ["0", "-2"])
def test_validate_args_rejects_non_positive_jobs(jobs):
    _expect_system_exit(["--jobs", jobs])


def test_validate_args_rejects_pvalue_without_null_distribution():
    _expect_system_exit(["--comparison-criterion", "p-value"])


def test_validate_args_rejects_out_of_range_lpd_for_sitega():
    _expect_system_exit(["-t", "sitega", "--lpd", "7"])


def test_validate_args_leaves_comparison_threshold_to_the_pipeline():
    args = _parse([])

    validate_args(create_arg_parser(), args)

    assert args.comparison_threshold is None


def test_check_dependencies_accepts_explicit_commands(tmp_path, monkeypatch):
    for env in (
        "HORDEMOTIFS_STREME_COMMAND",
        "HORDEMOTIFS_MEME_COMMAND",
        "HORDEMOTIFS_BAMM_COMMAND",
        "HORDEMOTIFS_SITEGA_COMMAND",
        "HORDEMOTIFS_DIMONT_JAR",
        "HORDEMOTIFS_SLIM_JAR",
    ):
        monkeypatch.delenv(env, raising=False)
    executable = tmp_path / "fake-tool"
    executable.touch()
    jar = tmp_path / "tool.jar"
    jar.touch()
    cases = {
        "streme": ["--streme-command", str(executable)],
        "meme": ["--meme-command", str(executable)],
        "bamm": [
            "--streme-command",
            str(executable),
            "--bamm-command",
            str(executable),
        ],
        "sitega": ["--sitega-command", str(executable)],
        "dimont": ["--java-command", str(executable), "--dimont-jar", str(jar)],
        "slim": ["--java-command", str(executable), "--slim-jar", str(jar)],
    }
    for tool, extra in cases.items():
        check_dependencies(_parse(["-t", tool, *extra]))


def test_check_dependencies_reports_missing_command():
    args = _parse(["-t", "streme", "--streme-command", "/nonexistent/streme"])

    with pytest.raises(RuntimeError, match="STREME dependency missing"):
        check_dependencies(args)


def test_check_dependencies_reports_missing_dimont_jar(tmp_path, monkeypatch):
    monkeypatch.delenv("HORDEMOTIFS_DIMONT_JAR", raising=False)
    executable = tmp_path / "java"
    executable.touch()
    args = _parse(["-t", "dimont", "--java-command", str(executable)])

    with pytest.raises(RuntimeError, match="HORDEMOTIFS_DIMONT_JAR"):
        check_dependencies(args)


def test_main_cli_help_exits_zero(capsys):
    with pytest.raises(SystemExit) as exc:
        main_cli(["--help"])

    assert exc.value.code == 0
    assert "usage:" in capsys.readouterr().out
