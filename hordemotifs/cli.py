"""
Command-line interface for HordeMotifs pipeline.
"""
import os
import sys
import shutil
import argparse
from typing import Dict, List, Any

from hordemotifs.pipeline import DeNovoPipeline
from hordemotifs.discovery import (
    StremeDiscoveryTool,
    BammDiscoveryTool,
    SitegaDiscoveryTool
)
from hordemotifs.evaluation import PerformanceEvaluator
from hordemotifs.comparison import (
    UniversalMotifComparator,
    TomtomComparator,
    MotaliComparator
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
    if '-' in s and ',' not in s:
        parts = s.split('-')
        if len(parts) != 3:
            raise ValueError(f"Invalid format for range: {s}. Expected 'start-end-step'")
        start, end, step = map(int, parts)
        return list(range(start, end + 1, step))
    elif ',' in s:
        return [int(x.strip()) for x in s.split(',')]
    else:
        return [int(s)]


def create_arg_parser() -> argparse.ArgumentParser:
    """Create and configure the command-line argument parser for HordeMotifs.
    
    Sets up all required and optional arguments for the motif discovery pipeline,
    including input files, discovery tool options, evaluation parameters, and
    motif comparison methods.
    
    Returns
    -------
    argparse.ArgumentParser
        Configured argument parser with all HordeMotifs options
    """
    parser = argparse.ArgumentParser(
        description='HordeMotifs: De novo motif discovery pipeline with odd/even bootstrap validation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
    Examples:
      # Basic PWM discovery with STREME
      hordemotifs peaks.fa background.fa promoters.fa output/ -t streme -l 8-20-4
    
      # Markov model-based motifs with different orders
      hordemotifs peaks.fa bg.fa promoters.fa output/ -t bamm -l 10-14-2 -o 1-4-1
    
      # SiteGA with custom LPD range
      hordemotifs peaks.fa bg.fa promoters.fa output/ -t sitega -l 10-16-2 --lpd 10-40-10
    
      # Single length value (when testing specific length)
      hordemotifs peaks.fa bg.fa promoters.fa output/ -t streme -l 12
    
      # Multiple specific values using comma-separated format
      hordemotifs peaks.fa bg.fa promoters.fa output/ -t bamm -l 10,12,14 -o 1,2,3
    
      # Custom comparison with Jaccard metric
      hordemotifs peaks.fa bg.fa promoters.fa output/ -c jaccard --c-metric jo --c-perm 5000
    
      # TomTom comparison with custom threshold
      hordemotifs peaks.fa bg.fa promoters.fa output/ -c tomtom --tomtom-pval 0.0001 --tomtom-metric pcc
    """
    )

    # ========== Required arguments ==========
    required = parser.add_argument_group('Required arguments')
    required.add_argument(
        'foreground',
        help='Path to the foreground FASTA file containing peak sequences'
    )
    required.add_argument(
        'background',
        help='Path to the background FASTA file'
    )
    required.add_argument(
        'promoters',
        help='Path to the promoter FASTA file used for threshold calculation'
    )
    required.add_argument(
        'output',
        help='Path to the output directory where results will be saved'
    )

    # ========== Discovery tool options ==========
    discovery = parser.add_argument_group('Motif discovery options')
    discovery.add_argument(
        '-t', '--tool',
        choices=['streme', 'bamm', 'sitega'],
        default='streme',
        help='De novo motif discovery tool to use (default: %(default)s)'
    )
    discovery.add_argument(
        '-n', '--nmotifs',
        type=int,
        default=5,
        help='Number of motifs to discover per run (default: %(default)s)'
    )
    discovery.add_argument(
        '-l', '--length',
        type=str,
        default='8-20-4',
        help='Range of motif lengths to discover. Format: \'start-end-step\' (e.g., 8-20-4), comma-separated list (e.g., 8,10,12), or a single value (default: %(default)s)'
    )
    discovery.add_argument(
        '-o', '--order',
        type=str,
        default='1-4-1',
        help='Range of Markov model orders. Format: \'start-end-step\' (e.g., 1-4-1), comma-separated list (e.g., 1,2,3), or a single value (default: %(default)s)'
    )
    discovery.add_argument(
        '--lpd',
        type=str,
        default='10-40-10',
        help='Range of locally positioned dinucleotide (LPD) distances for SiteGA. Format: \'start-end-step\', comma-separated list, or a single value (default: %(default)s)'
    )

    # ========== Evaluation options ==========
    evaluation = parser.add_argument_group('Evaluation options')
    evaluation.add_argument(
        '-f', '--fpr',
        type=float,
        default=0.001,
        help='False Positive Rate (FPR) threshold for partial AUC calculation (default: %(default)s)'
    )
    evaluation.add_argument(
        '-b', '--background-type',
        choices=['sites', 'peaks'],
        default='peaks',
        help='Method for background scoring. \'peaks\' uses the best site per sequence; \'sites\' uses all sites (default: %(default)s)'
    )
    evaluation.add_argument(
        '-m', '--metric',
        choices=['auROC', 'auPRC', 'pauROC', 'pauPRC'],
        default='pauROC',
        help='Performance metric used to select the best motifs (default: %(default)s)'
    )

    # ========== Comparison method options ==========
    comparison = parser.add_argument_group('Motif comparison options')
    comparison.add_argument(
        '-c', '--comparator',
        choices=['tomtom', 'continuous', 'motali'],
        default='tomtom',
        help='Method used for comparing discovered motifs (default: %(default)s)'
    )

    # TomTom options
    tomtom = parser.add_argument_group('TomTom comparator options')
    tomtom.add_argument(
        '--tomtom-metric',
        choices=['pcc', 'ed'],
        default='pcc',
        help='Distance metric for TomTom comparison. \'pcc\' is Pearson Correlation Coefficient; \'ed\' is Euclidean Distance (default: %(default)s)'
    )
    tomtom.add_argument(
        '--tomtom-pval',
        type=float,
        default=0.001,
        help='P-value threshold for TomTom Monte Carlo significance (default: %(default)s)'
    )
    tomtom.add_argument(
        '--tomtom-perm',
        type=int,
        default=1000,
        help='Number of permutations for TomTom significance testing (default: %(default)s)'
    )
    tomtom.add_argument(
        '--tomtom-permute-rows',
        action='store_true',
        help='Permute rows instead of columns when generating the null distribution'
    )
    tomtom.add_argument(
        '--tomtom-jobs',
        type=int,
        default=1,
        help='Number of parallel jobs for TomTom calculations (default: %(default)s)'
    )
    tomtom.add_argument(
        '--pfm-mode',
        action='store_true',
        help='If set, a Position Frequency Matrix (PFM) is derived for the model motifs by scanning sequences and constructing the PFM based on the top 5`%` of predicted binding sites'
    )

    continuous = parser.add_argument_group('Continuous comparator options')
    continuous.add_argument(
        '--c-metric',
        choices=['сj', 'co', 'corr'],
        default='сj',
        help='Metric for comparing motif models. Choices: cj (Continuous Jaccard), co (Continuous Overlap), corr (Pearson Correlation). (default: %(default)s)'
    )
    continuous.add_argument(
        '--c-perm',
        type=int,
        default=1000,
        help='Number of permutations for significance testing (default: %(default)s)'
    )
    continuous.add_argument(
        '--c-distortion',
        type=float,
        default=0.4,
        help='Distortion level applied during surrogate motif generation (default: %(default)s)'
    )
    continuous.add_argument(
        '--c-filter',
        choices=['score', 'p-value', 'none'],
        default='p-value',
        help='Criterion for filtering comparison results (default: %(default)s)'
    )
    continuous.add_argument(
        '--c-threshold',
        type=float,
        default=0.05,
        help='Numerical threshold for the filtering criterion (default: %(default)s)'
    )
    continuous.add_argument(
        '--c-search-range',
        type=int,
        default=10,
        help='Range to search for optimal offset alignment (default: %(default)s)'
    )
    continuous.add_argument(
        '--min-kernel-size',
        type=int,
        default=3,
        help='Minimum kernel size for convolution during surrogate generation. Used for `cj` and `co` options. (default: %(default)s)'
    )
    continuous.add_argument(
        '--max-kernel-size',
        type=int,
        default=11,
        help='Maximum kernel size for convolution during surrogate generation. Used for `cj` and `co` options. (default: %(default)s)'
    )
    continuous.add_argument(
        '--c-jobs',
        type=int,
        default=-1,
        help='Number of parallel jobs for comparisons. Used for `cj` and `co` options. Set to -1 to use all available cores (default: %(default)s)'
    )

    # Motali options
    motali = parser.add_argument_group('Motali comparator options')
    motali.add_argument(
        '--motali-threshold',
        type=float,
        default=0.95,
        help='Similarity threshold for Motali comparison (default: %(default)s)'
    )

    # ========== Other options ==========
    other = parser.add_argument_group('Other options')
    other.add_argument(
        '--seed',
        type=int,
        default=None,
        help='Random seed for reproducible results (default: %(default)s)'
    )
    other.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output logging'
    )

    return parser


def check_dependencies(tool: str) -> None:
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
    if tool == 'streme':
        if shutil.which("streme") is None:
            print("ERROR: STREME dependency missing. Please install MEME Suite from https://meme-suite.org/")
            sys.exit(1)

    elif tool == 'bamm':
        if shutil.which("BaMMmotif") is None:
            print("ERROR: BaMMmotif dependency missing. Please install BaMM from https://github.com/soedinglab/BaMMmotif")
            sys.exit(1)

    # Check Java for tools that need it
    if tool in ['chipmunk']:
        if shutil.which("java") is None:
            print("ERROR: Java dependency missing. Please install Java from https://www.java.com/")
            sys.exit(1)


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
    if args.tool == 'streme':
        return StremeDiscoveryTool(nmotifs=args.nmotifs)
    elif args.tool == 'bamm':
        return BammDiscoveryTool()
    elif args.tool == 'sitega':
        return SitegaDiscoveryTool(nmotifs=args.nmotifs)
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
    return PerformanceEvaluator(
        background_type=args.background_type
    )


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
        Instance of the appropriate comparator class (TomtomComparator,
        MotifComparator, or MotaliComparator)
    
    Raises
    ------
    ValueError
        If an unknown comparator name is specified in the arguments
    """
    if args.comparator == 'tomtom':
        return TomtomComparator(
            metric=args.tomtom_metric,
            n_permutations=args.tomtom_perm,
            permute_rows=args.tomtom_permute_rows,
            n_jobs=args.tomtom_jobs,
            seed=args.seed
        )

    elif args.comparator == 'continuous':
        filter_type = None if args.jo_filter == 'none' else args.jo_filter
        return UniversalMotifComparator(
            name=f"{args.comparator}_comparator",
            metric=args.c_metric,
            n_permutations=args.c_perm,
            distortion_level=args.c_distortion,
            n_jobs=args.c_jobs,
            seed=args.seed,
            filter_type=filter_type,
            filter_threshold=args.c_threshold,
            min_kernel_size = args.c_min_kernel_size,
            max_kernel_size = args.c_max_kernel_size,
            search_range = args.c_search_range

        )

    elif args.comparator == 'motali':
        return MotaliComparator(
            fasta_path=args.foreground,
            threshold=args.motali_threshold,
            tmp_directory=args.output
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
        params['length'] = parse_range(args.length)
    except ValueError as _:
        print(f"ERROR: Invalid format for --length: {args.length}. Expected 'start-end-step' or comma-separated values")
        sys.exit(1)

    # Parse order parameter (for Markov-based models)
    if args.tool in ['bamm']:
        try:
            params['order'] = parse_range(args.order)
        except ValueError as _:
            print(f"ERROR: Invalid format for --order: {args.order}. Expected 'start-end-step' or comma-separated values")
            sys.exit(1)

    # Parse LPD parameter (for SiteGA)
    if args.tool == 'sitega':
        try:
            params['lpd'] = parse_range(args.lpd)
        except ValueError as _:
            print(f"ERROR: Invalid format for --lpd: {args.lpd}. Expected 'start-end-step' or comma-separated values")
            sys.exit(1)

    if args.verbose:
        print("Discovery parameters:")
        for key, values in params.items():
            print(f"  {key}: {values}")

    return params


def main_cli():
    """Main command-line interface entry point for the HordeMotifs pipeline.
    
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
    check_dependencies(args.tool)

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    # Setup pipeline components
    if args.verbose:
        print("=" * 60)
        print("HordeMotifs De Novo Pipeline")
        print("=" * 60)
        print(f"Discovery tool: {args.tool}")
        print(f"Comparator: {args.comparator}")
        print(f"Metric: {args.metric}")
        print(f"FPR threshold: {args.fpr}")
        print("=" * 60)

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
    )

    pipeline.run(
        foreground_path=args.foreground,
        background_path=args.background,
        promoters_path=args.promoters,
        output_dir=args.output,
        discovery_params=discovery_params,
        metric=args.metric
    )

    if args.verbose:
        print("=" * 60)
        print("Pipeline completed successfully!")
        print(f"Results saved to: {args.output}")
        print("=" * 60)


if __name__ == '__main__':
    main_cli()
