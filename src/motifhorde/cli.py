"""Command-line interface for MotifHORDE."""

from __future__ import annotations

import argparse
import os
import shutil
from collections.abc import Callable
from dataclasses import dataclass

from .discovery import (
    BammDiscoveryTool,
    DimontDiscoveryTool,
    MemeDiscoveryTool,
    MotifDiscoveryTool,
    SitegaDiscoveryTool,
    SlimDiscoveryTool,
    StremeDiscoveryTool,
)
from .evaluation import PerformanceEvaluator
from .external import (
    DEFAULT_BAMM_COMMAND,
    DEFAULT_MEME_COMMAND,
    DEFAULT_SITEGA_COMMAND,
    DEFAULT_STREME_COMMAND,
    resolve_command,
    resolve_existing_path,
)
from .pipeline import DeNovoPipeline


def parse_range(value: str) -> list[int]:
    if "," in value:
        result = [int(item.strip()) for item in value.split(",")]
    elif value.count("-") == 2:
        start, end, step = (int(item) for item in value.split("-"))
        if step == 0:
            raise ValueError("range step must not be zero")
        result = list(range(start, end + (1 if step > 0 else -1), step))
    else:
        result = [int(value)]
    if not result:
        raise ValueError(f"empty range: {value}")
    return result


def create_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="MotifHORDE de novo motif discovery with odd/even bootstrap validation"
    )
    parser.add_argument("foreground")
    parser.add_argument("background")
    parser.add_argument("promoters")
    parser.add_argument("output")

    discovery = parser.add_argument_group("discovery")
    discovery.add_argument(
        "-t",
        "--tool",
        choices=["streme", "meme", "bamm", "dimont", "slim", "sitega"],
        default="streme",
    )
    discovery.add_argument("-n", "--nmotifs", type=int, default=5)
    discovery.add_argument("-l", "--length", default="8-20-4")
    discovery.add_argument("-o", "--order", default="1-4-1")
    discovery.add_argument(
        "--lpd",
        type=int,
        default=6,
        help="SiteGA maximum LPD distance (1..6)",
    )
    discovery.add_argument(
        "--features",
        default="10-40-10",
        help="SiteGA feature-count grid",
    )
    discovery.add_argument("--meme-command")
    discovery.add_argument("--streme-command")
    discovery.add_argument("--bamm-command")
    discovery.add_argument("--sitega-command")
    discovery.add_argument("--dimont-jar")
    discovery.add_argument("--slim-jar")
    discovery.add_argument("--java-command", default="java")
    discovery.add_argument("--java-xmx", default="4G")
    discovery.add_argument("--meme-objfun", default="classic")
    discovery.add_argument("--meme-mod", default="zoops")
    discovery.add_argument("--meme-minsites", type=int)
    discovery.add_argument("--meme-maxsites", type=int)
    discovery.add_argument("--meme-seed", type=int)
    discovery.add_argument("--jstacs-position-tag", default="position")
    discovery.add_argument("--jstacs-value-tag", default="value")
    discovery.add_argument("--jstacs-bg-order", type=int, default=-1)
    discovery.add_argument("--dimont-motif-order", type=int, default=0)
    discovery.add_argument("--dimont-ess", type=float, default=4.0)
    discovery.add_argument("--dimont-starts", type=int, default=20)
    discovery.add_argument("--slim-motif-order", type=int, default=-5)
    discovery.add_argument("--slim-starts", type=int, default=20)
    discovery.add_argument(
        "--slim-modify", action=argparse.BooleanOptionalAction, default=True
    )

    evaluation = parser.add_argument_group("evaluation")
    evaluation.add_argument("-f", "--fpr", type=float, default=0.001)
    evaluation.add_argument(
        "-b", "--background-type", choices=["sites", "peaks"], default="peaks"
    )
    evaluation.add_argument(
        "-m",
        "--metric",
        choices=["auROC", "auPRC", "pauROC", "pauPRC"],
        default="pauROC",
    )

    comparison = parser.add_argument_group("Mimosa profile comparison")
    comparison.add_argument(
        "--comparison-criterion", choices=["score", "p-value"], default="score"
    )
    comparison.add_argument("--comparison-threshold", type=float)
    comparison.add_argument("--mimosa-search-range", type=int, default=10)
    comparison.add_argument("--mimosa-null-distribution")

    other = parser.add_argument_group("other")
    other.add_argument("--jobs", type=int, default=1)
    other.add_argument("--seed", type=int)
    other.add_argument("-v", "--verbose", action="store_true")
    return parser


def _resolve_jobs(jobs: int) -> int:
    return (os.cpu_count() or 1) if jobs == -1 else jobs


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    if args.jobs == 0 or args.jobs < -1:
        parser.error("--jobs must be -1 or a positive integer")
    if args.nmotifs < 1:
        parser.error("--nmotifs must be positive")
    if not 0 <= args.fpr <= 1:
        parser.error("--fpr must lie in [0, 1]")
    if args.comparison_criterion == "p-value" and not args.mimosa_null_distribution:
        parser.error(
            "--comparison-criterion p-value requires --mimosa-null-distribution"
        )
    if args.tool == "sitega" and not 1 <= args.lpd <= 6:
        parser.error("--lpd must be between 1 and 6")


def _command_exists(command: str) -> bool:
    return os.path.isfile(command) or shutil.which(command) is not None


@dataclass(frozen=True)
class CommandRequirement:
    """An executable a discovery tool resolves from CLI/env/PATH."""

    label: str
    default: str
    cli_attr: str
    env_var: str | None


@dataclass(frozen=True)
class JarRequirement:
    """A JAR file a discovery tool resolves from CLI or an env variable."""

    label: str
    cli_attr: str
    env_var: str


TOOL_COMMANDS: dict[str, tuple[CommandRequirement, ...]] = {
    "streme": (
        CommandRequirement(
            "STREME",
            DEFAULT_STREME_COMMAND,
            "streme_command",
            "HORDEMOTIFS_STREME_COMMAND",
        ),
    ),
    "meme": (
        CommandRequirement(
            "MEME", DEFAULT_MEME_COMMAND, "meme_command", "HORDEMOTIFS_MEME_COMMAND"
        ),
    ),
    "bamm": (
        CommandRequirement(
            "STREME",
            DEFAULT_STREME_COMMAND,
            "streme_command",
            "HORDEMOTIFS_STREME_COMMAND",
        ),
        CommandRequirement(
            "BaMMmotif",
            DEFAULT_BAMM_COMMAND,
            "bamm_command",
            "HORDEMOTIFS_BAMM_COMMAND",
        ),
    ),
    "sitega": (
        CommandRequirement(
            "SiteGA",
            DEFAULT_SITEGA_COMMAND,
            "sitega_command",
            "HORDEMOTIFS_SITEGA_COMMAND",
        ),
    ),
    "dimont": (CommandRequirement("Java", "java", "java_command", None),),
    "slim": (CommandRequirement("Java", "java", "java_command", None),),
}

TOOL_JARS: dict[str, tuple[JarRequirement, ...]] = {
    "dimont": (JarRequirement("Dimont jar", "dimont_jar", "HORDEMOTIFS_DIMONT_JAR"),),
    "slim": (JarRequirement("SlimDimont jar", "slim_jar", "HORDEMOTIFS_SLIM_JAR"),),
}


def check_dependencies(args: argparse.Namespace) -> None:
    for requirement in TOOL_COMMANDS.get(args.tool, ()):
        command = resolve_command(
            requirement.default,
            explicit=getattr(args, requirement.cli_attr),
            env_var=requirement.env_var,
        )
        if not _command_exists(command):
            raise RuntimeError(f"{requirement.label} dependency missing: {command}")
    for requirement in TOOL_JARS.get(args.tool, ()):
        try:
            resolve_existing_path(
                getattr(args, requirement.cli_attr),
                requirement.env_var,
                requirement.label,
            )
        except FileNotFoundError as exc:
            raise RuntimeError(str(exc)) from exc


def _build_streme(args: argparse.Namespace) -> MotifDiscoveryTool:
    return StremeDiscoveryTool(args.streme_command)


def _build_meme(args: argparse.Namespace) -> MotifDiscoveryTool:
    return MemeDiscoveryTool(
        command=args.meme_command,
        objfun=args.meme_objfun,
        model=args.meme_mod,
        minsites=args.meme_minsites,
        maxsites=args.meme_maxsites,
        seed=args.meme_seed,
        threads=1,
    )


def _build_bamm(args: argparse.Namespace) -> MotifDiscoveryTool:
    return BammDiscoveryTool(args.bamm_command, args.streme_command)


def _build_dimont(args: argparse.Namespace) -> MotifDiscoveryTool:
    return DimontDiscoveryTool(
        jar_path=args.dimont_jar,
        java_command=args.java_command,
        java_xmx=args.java_xmx,
        threads=1,
        position_tag=args.jstacs_position_tag,
        value_tag=args.jstacs_value_tag,
        bg_order=args.jstacs_bg_order,
        motif_order=args.dimont_motif_order,
        ess=args.dimont_ess,
        starts=args.dimont_starts,
    )


def _build_slim(args: argparse.Namespace) -> MotifDiscoveryTool:
    return SlimDiscoveryTool(
        jar_path=args.slim_jar,
        java_command=args.java_command,
        java_xmx=args.java_xmx,
        threads=1,
        position_tag=args.jstacs_position_tag,
        value_tag=args.jstacs_value_tag,
        bg_order=args.jstacs_bg_order,
        motif_order=args.slim_motif_order,
        modify=args.slim_modify,
        starts=args.slim_starts,
    )


def _build_sitega(args: argparse.Namespace) -> MotifDiscoveryTool:
    return SitegaDiscoveryTool(args.sitega_command, 1, args.seed)


TOOL_BUILDERS: dict[str, Callable[[argparse.Namespace], MotifDiscoveryTool]] = {
    "streme": _build_streme,
    "meme": _build_meme,
    "bamm": _build_bamm,
    "dimont": _build_dimont,
    "slim": _build_slim,
    "sitega": _build_sitega,
}


def setup_discovery_tool(args: argparse.Namespace) -> MotifDiscoveryTool:
    try:
        builder = TOOL_BUILDERS[args.tool]
    except KeyError:
        raise ValueError(f"unknown discovery tool: {args.tool}") from None
    return builder(args)


def setup_discovery_params(args: argparse.Namespace) -> dict[str, list[int]]:
    params = {"length": parse_range(args.length)}
    if args.tool == "bamm":
        params["order"] = parse_range(args.order)
    if args.tool == "sitega":
        params["lpd"] = [args.lpd]
        params["features"] = parse_range(args.features)
    return params


def setup_evaluator(args: argparse.Namespace) -> PerformanceEvaluator:
    return PerformanceEvaluator(args.background_type)


def main_cli(argv: list[str] | None = None) -> None:
    parser = create_arg_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)
    check_dependencies(args)
    discovery_params = setup_discovery_params(args)
    if args.verbose:
        print(f"Discovery tool: {args.tool}")
        print(f"Discovery parameters: {discovery_params}")
        print(f"Jobs: {_resolve_jobs(args.jobs)}")

    pipeline = DeNovoPipeline(
        discovery_tool=setup_discovery_tool(args),
        evaluator=setup_evaluator(args),
        comparison_criterion=args.comparison_criterion,
        comparison_threshold=args.comparison_threshold,
        null_distribution=args.mimosa_null_distribution,
        fpr_threshold=args.fpr,
        number_of_motifs=args.nmotifs,
        jobs=_resolve_jobs(args.jobs),
        seed=args.seed,
        mimosa_search_range=args.mimosa_search_range,
    )
    pipeline.run(
        args.foreground,
        args.background,
        args.promoters,
        args.output,
        discovery_params,
        args.metric,
    )


if __name__ == "__main__":
    main_cli()
