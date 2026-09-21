from __future__ import annotations

from motifhorde.cli import create_arg_parser, setup_discovery_params


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
