module MotifHORDE

using ArgParse
using JSON3
using Logging
using Mimosa
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
