"""Command-line contract for the Julia runtime."""

function parse_range(value::AbstractString)
    text = strip(value)
    isempty(text) && throw(ArgumentError("range cannot be empty"))
    if occursin(',', text)
        return [parse(Int, strip(part)) for part in split(text, ',')]
    end
    matched = Base.match(r"^(-?\d+)-(-?\d+)-(-?\d+)$", text)
    if matched !== nothing
        start = parse(Int, something(matched.captures[1]))
        stop = parse(Int, something(matched.captures[2]))
        step = parse(Int, something(matched.captures[3]))
        step == 0 && throw(ArgumentError("range step cannot be zero"))
        return collect(start:step:stop)
    end
    return [parse(Int, text)]
end

function _cli_settings()
    settings = ArgParseSettings(;
        description="MotifHORDE: odd/even bootstrap motif discovery"
    )
    @add_arg_table! settings begin
        "foreground"
        help = "foreground FASTA"
        "background"
        help = "background FASTA"
        "promoters"
        help = "promoter FASTA"
        "output"
        help = "output directory"
        "-t", "--tool"
        arg_type = String
        default = "streme"
        range_tester =
            value -> value in ("streme", "meme", "bamm", "dimont", "slim", "sitega")
        "-n", "--nmotifs"
        arg_type = Int
        default = 5
        "-l", "--length"
        arg_type = String
        default = "8-20-4"
        "-o", "--order"
        arg_type = String
        default = "1-4-1"
        "--lpd"
        arg_type = String
        default = "10-40-10"
        "--streme-command"
        arg_type = String
        default = nothing
        "--meme-command"
        arg_type = String
        default = nothing
        "--bamm-command"
        arg_type = String
        default = nothing
        "--sitega-command"
        arg_type = String
        default = nothing
        "--dimont-jar"
        arg_type = String
        default = nothing
        "--slim-jar"
        arg_type = String
        default = nothing
        "--java-command"
        arg_type = String
        default = "java"
        "--java-xmx"
        arg_type = String
        default = "4G"
        "--fpr"
        arg_type = Float64
        default = 0.001
        "--background-type"
        arg_type = String
        default = "peaks"
        range_tester = value -> value in ("sites", "peaks")
        "-m", "--metric"
        arg_type = String
        default = "pauROC"
        range_tester = value -> value in VALIDATION_METRICS
        "-c", "--comparator"
        arg_type = String
        default = "tomtom"
        range_tester = value -> value in ("tomtom", "mimosa")
        "--tomtom-metric"
        arg_type = String
        default = "pcc"
        "--mimosa-metric"
        arg_type = String
        default = "co"
        "--comparison-criterion"
        arg_type = String
        default = "score"
        range_tester = value -> value in ("score", "p-value")
        "--comparison-threshold"
        arg_type = Float64
        default = nothing
        "--mimosa-search-range"
        arg_type = Int
        default = 10
        "--mimosa-null-distribution"
        arg_type = String
        default = nothing
        "--jobs"
        arg_type = Int
        default = 1
        "--seed"
        arg_type = Int
        default = nothing
        "-v", "--verbose"
        action = :store_true
    end
    return settings
end

function _setup_tool(args)
    if args["tool"] == "streme"
        return StremeDiscoveryTool(; command=args["streme-command"])
    elseif args["tool"] == "meme"
        return MemeDiscoveryTool(; command=args["meme-command"], threads=1)
    elseif args["tool"] == "bamm"
        return BammDiscoveryTool(;
            bamm_command=args["bamm-command"], streme_command=args["streme-command"]
        )
    elseif args["tool"] == "dimont"
        return DimontDiscoveryTool(;
            jar_path=args["dimont-jar"],
            java_command=args["java-command"],
            java_xmx=args["java-xmx"],
            threads=1,
        )
    elseif args["tool"] == "slim"
        return SlimDiscoveryTool(;
            jar_path=args["slim-jar"],
            java_command=args["java-command"],
            java_xmx=args["java-xmx"],
            threads=1,
        )
    elseif args["tool"] == "sitega"
        return SitegaDiscoveryTool(;
            executable=args["sitega-command"], threads=1, seed=args["seed"]
        )
    end
    return throw(ArgumentError("Unknown discovery tool: $(args["tool"])"))
end

function _setup_params(args)
    params = Dict{String,Any}("length" => parse_range(args["length"]))
    args["tool"] == "bamm" && (params["order"] = parse_range(args["order"]))
    args["tool"] == "sitega" && (params["lpd"] = parse_range(args["lpd"]))
    return params
end

function _validate_cli!(args)
    args["tool"] in ("streme", "meme", "bamm", "dimont", "slim", "sitega") ||
        throw(ArgumentError("Unsupported discovery tool: $(args["tool"])"))
    args["background-type"] in ("sites", "peaks") ||
        throw(ArgumentError("--background-type must be sites or peaks"))
    args["comparator"] in ("tomtom", "mimosa") ||
        throw(ArgumentError("--comparator must be tomtom or mimosa"))
    args["metric"] in VALIDATION_METRICS ||
        throw(ArgumentError("Unsupported validation metric: $(args["metric"])"))
    args["nmotifs"] > 0 || throw(ArgumentError("--nmotifs must be positive"))
    args["jobs"] == -1 ||
        args["jobs"] > 0 ||
        throw(ArgumentError("--jobs must be -1 or a positive integer"))
    all(>(0), _setup_params(args)["length"]) ||
        throw(ArgumentError("motif lengths must be positive"))
    return args
end

function configure_logging(verbose::Bool)
    level = verbose ? Logging.Info : Logging.Warn
    return global_logger(ConsoleLogger(stdout, level))
end
