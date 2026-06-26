"""
Command-line interface for MotifHORDE pipeline.
"""

import os
import sys
import shutil
import argparse
import logging
from typing import Dict, List, Any

from motifhorde.external import (
    DEFAULT_BAMM_COMMAND,
    DEFAULT_DIMONT_JAR,
    DEFAULT_MEME_COMMAND,
    DEFAULT_SLIM_JAR,
    DEFAULT_STREME_COMMAND,
    resolve_command,
    resolve_existing_path,
)
from motifhorde.pipeline import DeNovoPipeline
from motifhorde.discovery import (
    BammDiscoveryTool,
    DimontDiscoveryTool,
    MemeDiscoveryTool,
    SitegaDiscoveryTool,
    SlimDiscoveryTool,
    StremeDiscoveryTool,
)
from motifhorde.evaluation import PerformanceEvaluator
from motifhorde.comparison import (
    UniversalMotifComparator,
    TomtomComparator,
    default_threshold_for_criterion,
)


def parse_range(s: str) -> List[int]:
    """Parse a range string into a list of integers for bioinformatics parameter ranges.

    This function parses range strings in two formats:
    - Step format: 'start-end-step' (e.g., '8-20-4' becomes [8, 12, 16, 20])
    - Comma-separated format: 'value1,value2,value3' (e.g., '8,10,12')

    Parameters
    ----------
    s : str
        Range string in format 'start-end-step' or comma-separated values

    Returns
    -------
    List[int]
        List of integers parsed from the range string

    Raises
    ------
    ValueError
        If the range string format is invalid
    """
    if "-" in s and "," not in s:
        parts = s.split("-")
        if len(parts) != 3:
            raise ValueError(
                f"Invalid format for range: {s}. Expected 'start-end-step'"
            )
        start, end, step = map(int, parts)
        return list(range(start, end + 1, step))
    elif "," in s:
        return [int(x.strip()) for x in s.split(",")]
    else:
        return [int(s)]


def create_arg_parser() -> argparse.ArgumentParser:
    """Create and configure the command-line argument parser for MotifHORDE.

    Sets up all required and optional arguments for the motif discovery pipeline,
    including input files, discovery tool options, evaluation parameters, and
    motif comparison methods.

    Returns
    -------
    argparse.ArgumentParser
        Configured argument parser with all MotifHORDE options
    """
    parser = argparse.ArgumentParser(
        description="MotifHORDE: De novo motif discovery pipeline with odd/even bootstrap validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
    Examples:
      # Basic PWM discovery with STREME
      motifhorde peaks.fa background.fa promoters.fa output/ -t streme -l 8-20-4
    
      # Markov model-based motifs with different orders
      motifhorde peaks.fa bg.fa promoters.fa output/ -t bamm -l 10-14-2 -o 1-4-1
    
      # SiteGA with custom LPD range
      motifhorde peaks.fa bg.fa promoters.fa output/ -t sitega -l 10-16-2 --lpd 10-40-10
    
      # Single length value (when testing specific length)
      motifhorde peaks.fa bg.fa promoters.fa output/ -t streme -l 12
    
      # Multiple specific values using comma-separated format
      motifhorde peaks.fa bg.fa promoters.fa output/ -t bamm -l 10,12,14 -o 1,2,3
    
      # MIMOSA profile comparison
      motifhorde peaks.fa bg.fa promoters.fa output/ -c mimosa --mimosa-metric co

      # MIMOSA adjusted p-value filtering
      motifhorde peaks.fa bg.fa promoters.fa output/ -c mimosa --comparison-criterion p-value --mimosa-null-distribution profile-null.joblib

      # Matrix comparison with Pearson correlation
      motifhorde peaks.fa bg.fa promoters.fa output/ -c tomtom --tomtom-metric pcc
    """,
    )

    # ========== Required arguments ==========
    required = parser.add_argument_group("Required arguments")
    required.add_argument(
        "foreground", help="Path to the foreground FASTA file containing peak sequences"
    )
    required.add_argument("background", help="Path to the background FASTA file")
    required.add_argument(
        "promoters",
        help="Path to the promoter FASTA file used for threshold calculation",
    )
    required.add_argument(
        "output", help="Path to the output directory where results will be saved"
    )

    # ========== Discovery tool options ==========
    discovery = parser.add_argument_group("Motif discovery options")
    discovery.add_argument(
        "-t",
        "--tool",
        choices=["streme", "meme", "bamm", "dimont", "slim", "sitega"],
        default="streme",
        help="De novo motif discovery tool to use (default: %(default)s)",
    )
    discovery.add_argument(
        "-n",
        "--nmotifs",
        type=int,
        default=5,
        help="Number of motifs to discover per run (default: %(default)s)",
    )
    discovery.add_argument(
        "-l",
        "--length",
        type=str,
        default="8-20-4",
        help="Range of motif lengths to discover. Format: 'start-end-step' (e.g., 8-20-4), comma-separated list (e.g., 8,10,12), or a single value (default: %(default)s)",
    )
    discovery.add_argument(
        "-o",
        "--order",
        type=str,
        default="1-4-1",
        help="Range of Markov model orders. Format: 'start-end-step' (e.g., 1-4-1), comma-separated list (e.g., 1,2,3), or a single value (default: %(default)s)",
    )
    discovery.add_argument(
        "--lpd",
        type=str,
        default="10-40-10",
        help="Range of locally positioned dinucleotide (LPD) distances for SiteGA. Format: 'start-end-step', comma-separated list, or a single value (default: %(default)s)",
    )
    discovery.add_argument("--meme-command", default=None, help="MEME executable path")
    discovery.add_argument(
        "--streme-command", default=None, help="STREME executable path"
    )
    discovery.add_argument(
        "--bamm-command", default=None, help="BaMMmotif executable path"
    )
    discovery.add_argument("--dimont-jar", default=None, help="Dimont jar path")
    discovery.add_argument("--slim-jar", default=None, help="SlimDimont jar path")
    discovery.add_argument(
        "--java-command", default="java", help="Java executable path or command"
    )
    discovery.add_argument(
        "--java-xmx", default="4G", help="Java heap size for Jstacs tools"
    )
    discovery.add_argument(
        "--meme-objfun", default="classic", help="MEME objective function"
    )
    discovery.add_argument(
        "--meme-mod", default="zoops", help="MEME distribution model"
    )
    discovery.add_argument(
        "--meme-minsites", type=int, default=None, help="MEME minimum number of sites"
    )
    discovery.add_argument(
        "--meme-maxsites", type=int, default=None, help="MEME maximum number of sites"
    )
    discovery.add_argument(
        "--meme-seed", type=int, default=None, help="MEME random seed"
    )
    discovery.add_argument(
        "--jstacs-position-tag",
        default="position",
        help="Jstacs FASTA position annotation tag",
    )
    discovery.add_argument(
        "--jstacs-value-tag", default="value", help="Jstacs FASTA signal annotation tag"
    )
    discovery.add_argument(
        "--jstacs-bg-order", type=int, default=-1, help="Jstacs background Markov order"
    )
    discovery.add_argument(
        "--dimont-motif-order", type=int, default=0, help="Dimont motif Markov order"
    )
    discovery.add_argument(
        "--dimont-ess", type=float, default=4.0, help="Dimont equivalent sample size"
    )
    discovery.add_argument(
        "--dimont-starts", type=int, default=20, help="Dimont optimization starts"
    )
    discovery.add_argument(
        "--slim-motif-order",
        type=int,
        default=-5,
        help="SlimDimont motif order or negative distance",
    )
    discovery.add_argument(
        "--slim-starts", type=int, default=20, help="SlimDimont optimization starts"
    )
    discovery.add_argument(
        "--slim-modify",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enable SlimDimont shift adjustment",
    )

    # ========== Evaluation options ==========
    evaluation = parser.add_argument_group("Evaluation options")
    evaluation.add_argument(
        "-f",
        "--fpr",
        type=float,
        default=0.001,
        help="False Positive Rate (FPR) threshold for partial AUC calculation (default: %(default)s)",
    )
    evaluation.add_argument(
        "-b",
        "--background-type",
        choices=["sites", "peaks"],
        default="peaks",
        help="Method for background scoring. 'peaks' uses the best site per sequence; 'sites' uses all sites (default: %(default)s)",
    )
    evaluation.add_argument(
        "-m",
        "--metric",
        choices=["auROC", "auPRC", "pauROC", "pauPRC"],
        default="pauROC",
        help="Performance metric used to select the best motifs (default: %(default)s)",
    )

    # ========== Comparison method options ==========
    comparison = parser.add_argument_group("Motif comparison options")
    comparison.add_argument(
        "-c",
        "--comparator",
        choices=["tomtom", "mimosa"],
        default="tomtom",
        help="Method used for comparing discovered motifs (default: %(default)s)",
    )

    # TomTom options
    tomtom = parser.add_argument_group("TomTom comparator options")
    tomtom.add_argument(
        "--tomtom-metric",
        choices=["pcc", "ed"],
        default="pcc",
        help="Distance metric for TomTom comparison. 'pcc' is Pearson Correlation Coefficient; 'ed' is Euclidean Distance (default: %(default)s)",
    )
    tomtom.add_argument(
        "--pfm-mode",
        action="store_true",
        help="If set, derive PFMs by scanning sequences and using the top 5%% of predicted binding sites",
    )

    mimosa = parser.add_argument_group("MIMOSA comparator options")
    mimosa.add_argument(
        "--mimosa-metric",
        choices=["co", "co_rowwise", "dice", "dice_rowwise", "cosine"],
        default="co",
        help="Metric for comparing motif score profiles (default: %(default)s)",
    )
    mimosa.add_argument(
        "--comparison-criterion",
        choices=["score", "p-value"],
        default="score",
        help="Criterion for filtering comparison results. CLI p-value uses MIMOSA adj.p-value (default: %(default)s)",
    )
    mimosa.add_argument(
        "--comparison-threshold",
        type=float,
        default=None,
        help="Numerical threshold for the comparison criterion. Defaults to 0.9 for score and 0.05 for p-value",
    )
    mimosa.add_argument(
        "--mimosa-search-range",
        type=int,
        default=10,
        help="Range to search for optimal offset alignment (default: %(default)s)",
    )
    mimosa.add_argument(
        "--mimosa-null-distribution",
        default=None,
        help="Prepared MIMOSA null distribution path required for p-value filtering",
    )

    # ========== Other options ==========
    other = parser.add_argument_group("Other options")
    other.add_argument(
        "--jobs",
        type=int,
        default=1,
        help=(
            "Shared worker count. During bootstrap, independent discovery runs "
            "use this many processes and discovery tools run single-threaded "
            "inside each process. During comparison, supported comparators use "
            "this many internal jobs. Use -1 for all available cores "
            "(default: %(default)s)"
        ),
    )
    other.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible results (default: %(default)s)",
    )
    other.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose output logging"
    )

    return parser


def _resolve_jobs(jobs: int) -> int:
    if jobs == -1:
        return os.cpu_count() or 1
    return jobs


def configure_logging(verbose: bool) -> None:
    level = logging.INFO if verbose else logging.WARNING
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)-5s [%(name)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stdout,
        force=True,
    )


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    if args.jobs == 0 or args.jobs < -1:
        parser.error("--jobs must be -1 or a positive integer")

    if args.comparison_threshold is None:
        args.comparison_threshold = default_threshold_for_criterion(
            args.comparison_criterion
        )

    if args.comparison_criterion == "p-value" and args.comparator != "mimosa":
        parser.error("--comparison-criterion p-value is only supported with -c mimosa")

    if args.comparison_criterion == "p-value" and args.mimosa_null_distribution is None:
        parser.error(
            "--comparison-criterion p-value requires --mimosa-null-distribution"
        )


def _command_exists(command: str) -> bool:
    return os.path.exists(command) or shutil.which(command) is not None


def _dependency_error(message: str) -> None:
    print(f"ERROR: {message}")
    sys.exit(1)


def check_dependencies(args) -> None:
    """Check if required external dependencies are available in the system PATH.

    Verifies the presence of external tools needed for motif discovery based on
    the selected tool. Exits the program with an error if dependencies are missing.

    Parameters
    ----------
    tool : str
        Name of the discovery tool to check dependencies for ('streme', 'bamm', etc.)

    Raises
    ------
    SystemExit
        If required dependencies are not found in the system PATH
    """
    if args.tool == "streme":
        command = args.streme_command or resolve_command(
            "streme", DEFAULT_STREME_COMMAND, "HORDEMOTIFS_STREME_COMMAND"
        )
        if not _command_exists(command):
            _dependency_error(f"STREME dependency missing: {command}")
    elif args.tool == "meme":
        command = args.meme_command or resolve_command(
            "meme", DEFAULT_MEME_COMMAND, "HORDEMOTIFS_MEME_COMMAND"
        )
        if not _command_exists(command):
            _dependency_error(f"MEME dependency missing: {command}")
    elif args.tool == "bamm":
        streme = args.streme_command or resolve_command(
            "streme", DEFAULT_STREME_COMMAND, "HORDEMOTIFS_STREME_COMMAND"
        )
        bamm = args.bamm_command or resolve_command(
            DEFAULT_BAMM_COMMAND, DEFAULT_BAMM_COMMAND, "HORDEMOTIFS_BAMM_COMMAND"
        )
        if not _command_exists(streme):
            _dependency_error(
                f"STREME dependency missing for BaMM initialization: {streme}"
            )
        if not _command_exists(bamm):
            _dependency_error(f"BaMMmotif dependency missing: {bamm}")
    elif args.tool == "dimont":
        java = resolve_command(args.java_command)
        if not _command_exists(java):
            _dependency_error(f"Java dependency missing: {java}")
        try:
            resolve_existing_path(
                args.dimont_jar,
                "HORDEMOTIFS_DIMONT_JAR",
                DEFAULT_DIMONT_JAR,
                "Dimont jar",
            )
        except FileNotFoundError as exc:
            _dependency_error(str(exc))
    elif args.tool == "slim":
        java = resolve_command(args.java_command)
        if not _command_exists(java):
            _dependency_error(f"Java dependency missing: {java}")
        try:
            resolve_existing_path(
                args.slim_jar,
                "HORDEMOTIFS_SLIM_JAR",
                DEFAULT_SLIM_JAR,
                "SlimDimont jar",
            )
        except FileNotFoundError as exc:
            _dependency_error(str(exc))


def setup_discovery_tool(args) -> Any:
    """Initialize and return the appropriate motif discovery tool instance.

    Creates an instance of the specified discovery tool class based on command-line
    arguments, configuring it with the appropriate parameters.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments containing tool configuration

    Returns
    -------
    Any
        Instance of the appropriate discovery tool class (StremeDiscoveryTool,
        BammDiscoveryTool, or SitegaDiscoveryTool)

    Raises
    ------
    ValueError
        If an unknown tool name is specified in the arguments
    """
    discovery_threads = 1

    if args.tool == "streme":
        return StremeDiscoveryTool(nmotifs=args.nmotifs, command=args.streme_command)
    elif args.tool == "meme":
        return MemeDiscoveryTool(
            command=args.meme_command,
            objfun=args.meme_objfun,
            model=args.meme_mod,
            minsites=args.meme_minsites,
            maxsites=args.meme_maxsites,
            seed=args.meme_seed,
            threads=discovery_threads,
        )
    elif args.tool == "bamm":
        return BammDiscoveryTool(
            bamm_command=args.bamm_command, streme_command=args.streme_command
        )
    elif args.tool == "dimont":
        return DimontDiscoveryTool(
            jar_path=args.dimont_jar,
            java_command=args.java_command,
            java_xmx=args.java_xmx,
            threads=discovery_threads,
            position_tag=args.jstacs_position_tag,
            value_tag=args.jstacs_value_tag,
            bg_order=args.jstacs_bg_order,
            motif_order=args.dimont_motif_order,
            ess=args.dimont_ess,
            starts=args.dimont_starts,
        )
    elif args.tool == "slim":
        return SlimDiscoveryTool(
            jar_path=args.slim_jar,
            java_command=args.java_command,
            java_xmx=args.java_xmx,
            threads=discovery_threads,
            position_tag=args.jstacs_position_tag,
            value_tag=args.jstacs_value_tag,
            bg_order=args.jstacs_bg_order,
            motif_order=args.slim_motif_order,
            modify=args.slim_modify,
            starts=args.slim_starts,
        )
    elif args.tool == "sitega":
        return SitegaDiscoveryTool(
            nmotifs=args.nmotifs,
            threads=discovery_threads,
            seed=args.seed,
        )
    else:
        raise ValueError(f"Unknown tool: {args.tool}")


def setup_evaluator(args) -> PerformanceEvaluator:
    """Initialize and return the performance evaluator instance.

    Creates an instance of the PerformanceEvaluator class configured with
    parameters from the command-line arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments containing evaluator configuration

    Returns
    -------
    PerformanceEvaluator
        Configured performance evaluator instance
    """
    return PerformanceEvaluator(background_type=args.background_type)


def setup_comparator(args) -> Any:
    """Initialize and return the appropriate motif comparison tool instance.

    Creates an instance of the specified comparator class based on command-line
    arguments, configuring it with the appropriate parameters for motif comparison.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments containing comparator configuration

    Returns
    -------
    Any
        Instance of the appropriate comparator class.

    Raises
    ------
    ValueError
        If an unknown comparator name is specified in the arguments
    """
    if args.comparator == "tomtom":
        return TomtomComparator(
            metric=args.tomtom_metric,
            n_jobs=args.jobs,
            seed=args.seed,
            pfm_mode=args.pfm_mode,
            comparison_criterion=args.comparison_criterion,
            comparison_threshold=args.comparison_threshold,
        )

    elif args.comparator == "mimosa":
        return UniversalMotifComparator(
            name="mimosa_comparator",
            metric=args.mimosa_metric,
            n_jobs=args.jobs,
            seed=args.seed,
            comparison_criterion=args.comparison_criterion,
            comparison_threshold=args.comparison_threshold,
            search_range=args.mimosa_search_range,
            pvalue=args.comparison_criterion == "p-value",
            null_distribution=args.mimosa_null_distribution,
        )

    else:
        raise ValueError(f"Unknown comparator: {args.comparator}")


def setup_discovery_params(args) -> Dict[str, List[Any]]:
    """Create a dictionary of discovery parameters based on tool and arguments.

    Parses command-line arguments to generate a parameter dictionary for the
    motif discovery process, handling different parameter types based on the
    selected discovery tool.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments containing discovery parameter specifications

    Returns
    -------
    Dict[str, List[Any]]
        Dictionary containing discovery parameters like:
        {'length': [8, 12, 16, 20], 'order': [1, 2, 3, 4]}
    """
    params = {}

    # Parse length parameter (common for all tools)
    try:
        params["length"] = parse_range(args.length)
    except ValueError as _:
        print(
            f"ERROR: Invalid format for --length: {args.length}. Expected 'start-end-step' or comma-separated values"
        )
        sys.exit(1)

    # Parse order parameter (for Markov-based models)
    if args.tool in ["bamm"]:
        try:
            params["order"] = parse_range(args.order)
        except ValueError as _:
            print(
                f"ERROR: Invalid format for --order: {args.order}. Expected 'start-end-step' or comma-separated values"
            )
            sys.exit(1)

    # Parse LPD parameter (for SiteGA)
    if args.tool == "sitega":
        try:
            params["lpd"] = parse_range(args.lpd)
        except ValueError as _:
            print(
                f"ERROR: Invalid format for --lpd: {args.lpd}. Expected 'start-end-step' or comma-separated values"
            )
            sys.exit(1)

    if args.verbose:
        params_text = " | ".join(
            f"{key}={','.join(str(value) for value in values)}"
            for key, values in sorted(params.items())
        )
        logging.getLogger("pipeline").info("params | %s", params_text)

    return params


def main_cli():
    """Main command-line interface entry point for the MotifHORDE pipeline.

    Orchestrates the complete motif discovery pipeline by parsing command-line
    arguments, validating inputs, setting up pipeline components, and executing
    the discovery process.
    """
    # Parse arguments
    parser = create_arg_parser()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    validate_args(parser, args)
    configure_logging(args.verbose)

    # Validate inputs
    if not os.path.exists(args.foreground):
        print(f"ERROR: Foreground FASTA file not found: {args.foreground}")
        sys.exit(1)

    if not os.path.exists(args.background):
        print(f"ERROR: Background FASTA file not found: {args.background}")
        sys.exit(1)

    if not os.path.exists(args.promoters):
        print(f"ERROR: Promoter FASTA file not found: {args.promoters}")
        sys.exit(1)

    # Check dependencies
    check_dependencies(args)

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    # Setup pipeline components
    resolved_jobs = _resolve_jobs(args.jobs)
    pipeline_logger = logging.getLogger("pipeline")
    pipeline_logger.info(
        "start | tool=%s | comparator=%s | metric=%s | fpr=%s | jobs=%s | seed=%s",
        args.tool,
        args.comparator,
        args.metric,
        args.fpr,
        resolved_jobs,
        args.seed,
    )
    pipeline_logger.info(
        "workers | bootstrap_discovery=%s | discovery_threads=1 | comparator_jobs=%s",
        resolved_jobs,
        args.jobs,
    )

    discovery_tool = setup_discovery_tool(args)
    evaluator = setup_evaluator(args)
    comparator = setup_comparator(args)
    discovery_params = setup_discovery_params(args)

    # Run pipeline
    pipeline = DeNovoPipeline(
        discovery_tool=discovery_tool,
        evaluator=evaluator,
        comparator=comparator,
        fpr_threshold=args.fpr,
        number_of_motifs=args.nmotifs,
        jobs=_resolve_jobs(args.jobs),
        seed=args.seed,
    )

    pipeline.run(
        foreground_path=args.foreground,
        background_path=args.background,
        promoters_path=args.promoters,
        output_dir=args.output,
        discovery_params=discovery_params,
        metric=args.metric,
    )


if __name__ == "__main__":
    main_cli()
