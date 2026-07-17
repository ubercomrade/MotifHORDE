module MotifHORDE

using ArgParse: @add_arg_table!, ArgParseSettings, parse_args
using JSON3
using Logging
using Mimosa: Mimosa
using Printf
using Random
using Statistics

include("models.jl")
include("fasta.jl")
include("external.jl")
include("discovery.jl")
include("comparison.jl")
include("evaluation.jl")
include("bootstrap.jl")
include("pipeline.jl")
include("cli.jl")

function (@main)(arguments::Vector{String}=copy(ARGS))
    parsed = parse_args(args, _cli_settings())
    _validate_cli!(parsed)
    criterion = parsed["comparison-criterion"]
    parsed["comparator"] == "mimosa" ||
        criterion == "score" ||
        throw(ArgumentError("p-value criterion requires --comparator mimosa"))
    criterion == "score" ||
        parsed["mimosa-null-distribution"] !== nothing ||
        throw(ArgumentError("p-value criterion requires --mimosa-null-distribution"))
    configure_logging(parsed["verbose"])
    tool = _setup_tool(parsed)
    evaluator = PerformanceEvaluator(; background_type=Symbol(parsed["background-type"]))
    threshold = if parsed["comparison-threshold"] === nothing
        default_threshold_for_criterion(criterion)
    else
        parsed["comparison-threshold"]
    end
    comparator = if parsed["comparator"] == "mimosa"
        UniversalMotifComparator(;
            metric=Symbol(parsed["mimosa-metric"]),
            comparison_criterion=Symbol(criterion),
            comparison_threshold=threshold,
            search_range=parsed["mimosa-search-range"],
            null_distribution=parsed["mimosa-null-distribution"],
        )
    else
        TomtomComparator(;
            metric=Symbol(parsed["tomtom-metric"]),
            comparison_criterion=Symbol(criterion),
            comparison_threshold=threshold,
        )
    end
    jobs = parsed["jobs"] == -1 ? Threads.nthreads() : parsed["jobs"]
    config = PipelineConfig(parsed["fpr"], parsed["nmotifs"], jobs, parsed["seed"])
    run_pipeline(
        tool,
        evaluator,
        comparator,
        parsed["foreground"],
        parsed["background"],
        parsed["promoters"],
        parsed["output"],
        _setup_params(parsed);
        metric=parsed["metric"],
        config=config,
    )
    return 0
end

export BootstrapTask,
    Comparison,
    DimontDiscoveryTool,
    MotifDiscoveryTool,
    PerformanceEvaluator,
    PipelineConfig,
    SitegaDiscoveryTool,
    SlimDiscoveryTool,
    StremeDiscoveryTool,
    MemeDiscoveryTool,
    BammDiscoveryTool,
    TomtomComparator,
    UniversalMotifComparator,
    all_valid_scores,
    best_scores,
    build_bootstrap_tasks,
    build_dimont_args,
    build_sitega_args,
    build_slim_args,
    build_streme_args,
    comparison_column_for_criterion,
    default_threshold_for_criterion,
    discover,
    evaluate,
    filter_similar,
    format_params,
    parse_range,
    parameter_grid,
    read_fasta,
    read_model,
    reconstruct_pfm,
    run_pipeline,
    run_checked,
    select_sequence_rows,
    write_fasta,
    write_jstacs_fasta,
    write_meme,
    write_model

end
