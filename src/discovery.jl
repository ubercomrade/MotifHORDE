"""External discovery adapters with one ordered model-return contract."""

abstract type MotifDiscoveryTool end

discover(::MotifDiscoveryTool, args...; kwargs...) = throw(MethodError(discover, args))

function _require_length(kwargs)
    haskey(kwargs, :length) ||
        throw(ArgumentError("Parameter 'length' is required for discovery."))
    length = Int(kwargs[:length])
    length > 0 || throw(ArgumentError("Parameter 'length' must be positive."))
    return length
end

function _expected_length(model, expected::Union{Nothing,Int})
    expected === nothing || motif_length(model) == expected || return false
    return true
end

function _read_indexed_models(
    path::AbstractString,
    format::Symbol,
    prefix::AbstractString,
    number_of_motifs::Int;
    expected_length::Union{Nothing,Int}=nothing,
)
    isfile(path) && filesize(path) > 0 || return Mimosa.AbstractMotifModel[]
    models = Mimosa.AbstractMotifModel[]
    for index in 0:(number_of_motifs - 1)
        model = try
            read_model(path, format; index=index)
        catch error
            error isa BoundsError || error isa ArgumentError || break
            break
        end
        _expected_length(model, expected_length) || continue
        push!(models, rename_model(model, "$prefix-$(length(models) + 1)"))
    end
    return models
end

function _meme_width(header::AbstractString)
    matched = Base.match(r"\bw\s*=\s*(\d+)", header)
    return matched === nothing ? 0 : parse(Int, something(matched.captures[1]))
end

function _meme_report_name(line::AbstractString)
    parts = split(strip(line))
    length(parts) >= 2 && parts[1] == "MOTIF" || return nothing
    return length(parts) >= 3 && startswith(parts[3], "MEME-") ? parts[3] : parts[2]
end

function _is_probability_row(line::AbstractString)
    values = split(strip(line))
    length(values) == 4 || return false
    return all(tryparse(Float64, value) !== nothing for value in values)
end

"""Turn MEME's verbose text report into the compact reader input."""
function normalize_meme_report(text::AbstractString)
    lines = split(text, '\n'; keepempty=true)
    records = Tuple{String,String,Vector{String}}[]
    current_name = nothing
    index = 1
    while index <= length(lines)
        name = _meme_report_name(lines[index])
        name !== nothing && (current_name = name)
        if startswith(lines[index], "letter-probability matrix:")
            width = _meme_width(lines[index])
            stop = index + width
            if current_name !== nothing && width > 0 && stop <= length(lines)
                rows = String.(lines[(index + 1):stop])
                if all(_is_probability_row, rows)
                    push!(
                        records, (String(current_name), strip(lines[index]), strip.(rows))
                    )
                    index = stop
                end
            end
        end
        index += 1
    end
    isempty(records) && return String(text)
    output = ["MEME version 4", "", "ALPHABET= ACGT", "", "strands: + -", ""]
    for (name, header, rows) in records
        append!(output, ["MOTIF $name", header])
        append!(output, rows)
        push!(output, "")
    end
    return join(output, '\n')
end

struct StremeDiscoveryTool <: MotifDiscoveryTool
    command::Union{Nothing,String}
end
function StremeDiscoveryTool(; command=nothing)
    return StremeDiscoveryTool(command === nothing ? nothing : String(command))
end

function build_streme_args(
    command, foreground, background, output_dir, length, number_of_motifs
)
    return String[
        command,
        "--p",
        foreground,
        "--n",
        background,
        "--objfun",
        "de",
        "--minw",
        string(length),
        "--maxw",
        string(length),
        "-nmotifs",
        string(number_of_motifs),
        "--text",
    ]
end

function discover(
    tool::StremeDiscoveryTool,
    foreground::AbstractString,
    background::AbstractString,
    output_dir::AbstractString,
    number_of_motifs::Integer;
    length=nothing,
    kwargs...,
)
    width = length === nothing ? _require_length((; kwargs...)) : Int(length)
    mkpath(output_dir)
    command = resolve_command(
        tool.command, DEFAULT_STREME_COMMAND; env_var="HORDEMOTIFS_STREME_COMMAND"
    )
    result = run_checked(
        build_streme_args(
            command, foreground, background, output_dir, width, Int(number_of_motifs)
        ),
    )
    path = joinpath(output_dir, "motifs.meme")
    open(path, "w") do io
        return write(io, result.stdout)
    end
    return _read_indexed_models(
        path, :meme, "Streme", Int(number_of_motifs); expected_length=width
    )
end

struct MemeDiscoveryTool <: MotifDiscoveryTool
    command::Union{Nothing,String}
    objfun::String
    model::String
    minsites::Union{Nothing,Int}
    maxsites::Union{Nothing,Int}
    seed::Union{Nothing,Int}
    threads::Union{Nothing,Int}
end

function MemeDiscoveryTool(;
    command=nothing,
    objfun="classic",
    model="zoops",
    minsites=nothing,
    maxsites=nothing,
    seed=nothing,
    threads=nothing,
)
    return MemeDiscoveryTool(
        command === nothing ? nothing : String(command),
        String(objfun),
        String(model),
        minsites,
        maxsites,
        seed,
        threads,
    )
end

function build_meme_args(
    tool::MemeDiscoveryTool, command, foreground, background, length, number_of_motifs
)
    args = String[
        command,
        foreground,
        "-dna",
        "-revcomp",
        "-nmotifs",
        string(number_of_motifs),
        "-minw",
        string(length),
        "-maxw",
        string(length),
        "-nomatrim",
        "-text",
        "-objfun",
        tool.objfun,
        "-mod",
        tool.model,
    ]
    tool.objfun == "classic" || append!(args, ["-neg", background])
    for (flag, value) in (
        ("-minsites", tool.minsites),
        ("-maxsites", tool.maxsites),
        ("-seed", tool.seed),
        ("-p", tool.threads),
    )
        value === nothing || append!(args, [flag, string(value)])
    end
    return args
end

function discover(
    tool::MemeDiscoveryTool,
    foreground::AbstractString,
    background::AbstractString,
    output_dir::AbstractString,
    number_of_motifs::Integer;
    length=nothing,
    kwargs...,
)
    width = length === nothing ? _require_length((; kwargs...)) : Int(length)
    mkpath(output_dir)
    command = resolve_command(
        tool.command, DEFAULT_MEME_COMMAND; env_var="HORDEMOTIFS_MEME_COMMAND"
    )
    result = run_checked(
        build_meme_args(tool, command, foreground, background, width, Int(number_of_motifs))
    )
    path = joinpath(output_dir, "motifs.meme")
    if !isempty(result.stdout)
        open(path, "w") do io
            return write(io, normalize_meme_report(result.stdout))
        end
    elseif isfile(joinpath(output_dir, "meme.txt"))
        path = joinpath(output_dir, "meme.txt")
    end
    return _read_indexed_models(
        path, :meme, "Meme", Int(number_of_motifs); expected_length=width
    )
end

struct BammDiscoveryTool <: MotifDiscoveryTool
    bamm_command::Union{Nothing,String}
    streme_command::Union{Nothing,String}
end
function BammDiscoveryTool(; bamm_command=nothing, streme_command=nothing)
    return BammDiscoveryTool(bamm_command, streme_command)
end

function build_bamm_args(command, output_dir, foreground, background, meme_path, order)
    return String[
        command,
        output_dir,
        foreground,
        "--PWMFile",
        meme_path,
        "--EM",
        "--order",
        string(order),
        "--Order",
        string(order),
        "--basename",
        "bamm",
        "--negSeqFile",
        background,
    ]
end

function discover(
    tool::BammDiscoveryTool,
    foreground::AbstractString,
    background::AbstractString,
    output_dir::AbstractString,
    number_of_motifs::Integer;
    length=nothing,
    order=2,
    kwargs...,
)
    width = length === nothing ? _require_length((; kwargs...)) : Int(length)
    mkpath(output_dir)
    streme = StremeDiscoveryTool(; command=tool.streme_command)
    discover(streme, foreground, background, output_dir, number_of_motifs; length=width)
    meme_path = joinpath(output_dir, "motifs.meme")
    isfile(meme_path) || return Mimosa.AbstractMotifModel[]
    bamm = resolve_command(
        tool.bamm_command, DEFAULT_BAMM_COMMAND; env_var="HORDEMOTIFS_BAMM_COMMAND"
    )
    run_checked(
        build_bamm_args(bamm, output_dir, foreground, background, meme_path, Int(order))
    )
    models = Mimosa.AbstractMotifModel[]
    for index in 1:Int(number_of_motifs)
        path = joinpath(output_dir, "bamm_motif_$(index).ihbcp")
        isfile(path) || continue
        model = read_model(path, :bamm; order=Int(order))
        _expected_length(model, width) || continue
        push!(models, rename_model(model, "Bamm-$(length(models) + 1)"))
    end
    return models
end

function _existing_xml_paths(output_dir::AbstractString, patterns::Vector{String})
    files = String[]
    for name in readdir(output_dir)
        any(occursin(Regex(replace(pattern, "*" => ".*")), name) for pattern in patterns) &&
            push!(files, joinpath(output_dir, name))
    end
    return sort(filter(isfile, files))
end

function build_dimont_args(
    java,
    xmx,
    jar,
    output_dir,
    data_name,
    length,
    position_tag,
    value_tag,
    bg_order,
    motif_order,
    ess,
    starts,
    threads,
)
    args = String[
        java,
        "-Djava.awt.headless=true",
        "-Xmx$xmx",
        "-jar",
        jar,
        "home=$output_dir",
        "data=$data_name",
        "infix=dimont",
        "position=$position_tag",
        "value=$value_tag",
        "motifWidth=$(length)",
        "motifOrder=$(motif_order)",
        "bgOrder=$(bg_order)",
        "ess=$(ess)",
        "starts=$(starts)",
    ]
    threads === nothing || push!(args, "threads=$(threads)")
    return args
end

function build_slim_args(
    java,
    xmx,
    jar,
    output_dir,
    data_name,
    length,
    position_tag,
    value_tag,
    bg_order,
    motif_order,
    modify,
    starts,
    threads,
)
    args = String[
        java,
        "-Djava.awt.headless=true",
        "-Xmx$xmx",
        "-jar",
        jar,
        "home=$output_dir",
        "data=$data_name",
        "infix=slim",
        "position=$position_tag",
        "value=$value_tag",
        "motifWidth=$(length)",
        "motifOrder=$(motif_order)",
        "bgOrder=$(bg_order)",
        "starts=$(starts)",
        "modify=$(modify ? "true" : "false")",
    ]
    threads === nothing || push!(args, "threads=$(threads)")
    return args
end

struct DimontDiscoveryTool <: MotifDiscoveryTool
    jar_path::Union{Nothing,String}
    java_command::String
    java_xmx::String
    threads::Union{Nothing,Int}
    position_tag::String
    value_tag::String
    bg_order::Int
    motif_order::Int
    ess::Float64
    starts::Int
end

function DimontDiscoveryTool(;
    jar_path=nothing,
    java_command="java",
    java_xmx="4G",
    threads=nothing,
    position_tag="position",
    value_tag="value",
    bg_order=-1,
    motif_order=0,
    ess=4.0,
    starts=20,
)
    return DimontDiscoveryTool(
        jar_path,
        String(java_command),
        String(java_xmx),
        threads,
        String(position_tag),
        String(value_tag),
        bg_order,
        motif_order,
        Float64(ess),
        starts,
    )
end

function _discover_jstacs(
    tool, format, prefix, foreground, output_dir, number_of_motifs, width
)
    mkpath(output_dir)
    jar = resolve_existing_path(
        tool.jar_path,
        format == :dimont ? "HORDEMOTIFS_DIMONT_JAR" : "HORDEMOTIFS_SLIM_JAR",
        joinpath(
            dirname(@__DIR__),
            "artifacts",
            format == :dimont ? "Dimont.jar" : "SlimDimont.jar",
        ),
        format == :dimont ? "Dimont jar" : "SlimDimont jar",
    )
    data_name = "train.annot.fa"
    write_jstacs_fasta(
        foreground,
        joinpath(output_dir, data_name);
        position_tag=tool.position_tag,
        value_tag=tool.value_tag,
    )
    command = resolve_command(tool.java_command, tool.java_command)
    args = if format == :dimont
        build_dimont_args(
            command,
            tool.java_xmx,
            jar,
            output_dir,
            data_name,
            width,
            tool.position_tag,
            tool.value_tag,
            tool.bg_order,
            tool.motif_order,
            tool.ess,
            tool.starts,
            tool.threads,
        )
    else
        build_slim_args(
            command,
            tool.java_xmx,
            jar,
            output_dir,
            data_name,
            width,
            tool.position_tag,
            tool.value_tag,
            tool.bg_order,
            tool.motif_order,
            tool.modify,
            tool.starts,
            tool.threads,
        )
    end
    run_checked(args)
    paths = _existing_xml_paths(
        output_dir,
        if format == :dimont
            [".*dimont.*\\.xml", ".*Dimont.*\\.xml", ".*\\.xml"]
        else
            [".*slim.*\\.xml", ".*Slim.*\\.xml", ".*\\.xml"]
        end,
    )
    models = Mimosa.AbstractMotifModel[]
    for path in paths
        length(models) >= number_of_motifs && break
        model = try
            read_model(path, format)
        catch error
            @warn "Skipping invalid Jstacs output" path exception=(error, catch_backtrace())
            continue
        end
        _expected_length(model, width) || continue
        push!(models, rename_model(model, "$prefix-$(length(models) + 1)"))
    end
    return models
end

struct SlimDiscoveryTool <: MotifDiscoveryTool
    jar_path::Union{Nothing,String}
    java_command::String
    java_xmx::String
    threads::Union{Nothing,Int}
    position_tag::String
    value_tag::String
    bg_order::Int
    motif_order::Int
    modify::Bool
    starts::Int
end

function SlimDiscoveryTool(;
    jar_path=nothing,
    java_command="java",
    java_xmx="4G",
    threads=nothing,
    position_tag="position",
    value_tag="value",
    bg_order=-1,
    motif_order=-5,
    modify=true,
    starts=20,
)
    return SlimDiscoveryTool(
        jar_path,
        String(java_command),
        String(java_xmx),
        threads,
        String(position_tag),
        String(value_tag),
        bg_order,
        motif_order,
        Bool(modify),
        starts,
    )
end

function discover(
    tool::DimontDiscoveryTool,
    foreground,
    background,
    output_dir,
    number_of_motifs;
    length=nothing,
    kwargs...,
)
    width = length === nothing ? _require_length((; kwargs...)) : Int(length)
    return _discover_jstacs(
        tool, :dimont, "Dimont", foreground, output_dir, Int(number_of_motifs), width
    )
end

function discover(
    tool::SlimDiscoveryTool,
    foreground,
    background,
    output_dir,
    number_of_motifs;
    length=nothing,
    kwargs...,
)
    width = length === nothing ? _require_length((; kwargs...)) : Int(length)
    return _discover_jstacs(
        tool, :slim, "Slim", foreground, output_dir, Int(number_of_motifs), width
    )
end

struct SitegaDiscoveryTool <: MotifDiscoveryTool
    executable::Union{Nothing,String}
    threads::Union{Nothing,Int}
    seed::Union{Nothing,Int}
end
function SitegaDiscoveryTool(; executable=nothing, threads=nothing, seed=nothing)
    return SitegaDiscoveryTool(executable, threads, seed)
end

function build_sitega_args(
    command,
    foreground,
    background,
    output_dir,
    manifest,
    length,
    lpd,
    motifs,
    seed,
    threads,
)
    args = String[
        command,
        "--foreground",
        foreground,
        "--background",
        background,
        "--output",
        output_dir,
        "--manifest",
        manifest,
        "--length",
        string(length),
        "--lpd",
        string(lpd),
        "--motifs",
        string(motifs),
    ]
    seed === nothing || append!(args, ["--seed", string(seed)])
    threads === nothing || append!(args, ["--threads", string(threads)])
    return args
end

function _manifest_value(manifest, key::AbstractString)
    haskey(manifest, key) || throw(ArgumentError("SiteGA manifest is missing '$key'."))
    return manifest[key]
end

function _read_sitega_manifest(path, output_dir, requested, width)
    manifest = try
        JSON3.read(read(path, String))
    catch error
        throw(ArgumentError("Invalid SiteGA manifest $path: $(sprint(showerror, error))"))
    end
    _manifest_value(manifest, "schema_version") == 1 ||
        throw(ArgumentError("Unsupported SiteGA manifest schema_version."))
    _manifest_value(manifest, "status") == "success" || throw(
        ArgumentError("SiteGA reported failure: $(_manifest_value(manifest, "message"))"),
    )
    entries = _manifest_value(manifest, "models")
    models = Mimosa.AbstractMotifModel[]
    for (index, entry) in enumerate(entries)
        length(models) >= requested && break
        model_path = String(_manifest_value(entry, "model_file"))
        candidate = normpath(joinpath(output_dir, model_path))
        startswith(
            candidate, normpath(output_dir) * string(Base.Filesystem.path_separator)
        ) || throw(ArgumentError("SiteGA manifest model path escapes output directory."))
        isfile(candidate) || throw(ArgumentError("SiteGA model file not found: $candidate"))
        model_type = Symbol(lowercase(String(_manifest_value(entry, "model_type"))))
        model = read_model(candidate, model_type)
        declared_length = Int(_manifest_value(entry, "length"))
        declared_length == motif_length(model) ||
            throw(ArgumentError("SiteGA manifest length disagrees with model file."))
        declared_length == width || continue
        name = String(get(entry, "name", "Sitega-$index"))
        push!(models, rename_model(model, name))
    end
    return models
end

function discover(
    tool::SitegaDiscoveryTool,
    foreground,
    background,
    output_dir,
    number_of_motifs;
    length=nothing,
    lpd=20,
    seed=tool.seed,
    kwargs...,
)
    width = length === nothing ? _require_length((; kwargs...)) : Int(length)
    mkpath(output_dir)
    command = resolve_command(
        tool.executable, DEFAULT_SITEGA_COMMAND; env_var="HORDEMOTIFS_SITEGA_COMMAND"
    )
    manifest = joinpath(output_dir, "sitega.manifest.json")
    args = build_sitega_args(
        command,
        foreground,
        background,
        output_dir,
        manifest,
        width,
        Int(lpd),
        Int(number_of_motifs),
        seed,
        tool.threads,
    )
    run_checked(args)
    isfile(manifest) || throw(ArgumentError("SiteGA did not create manifest: $manifest"))
    return _read_sitega_manifest(manifest, output_dir, Int(number_of_motifs), width)
end
