"""Safe resolution and execution of external programs."""

const DEFAULT_MEME_COMMAND = "meme"
const DEFAULT_STREME_COMMAND = "streme"
const DEFAULT_BAMM_COMMAND = "BaMMmotif"

function resolve_command(
    cli_value::Union{Nothing,AbstractString},
    command_name::AbstractString;
    env_var::Union{Nothing,AbstractString}=nothing,
    fallback::Union{Nothing,AbstractString}=nothing,
)
    cli_value !== nothing && !isempty(cli_value) && return String(cli_value)
    if env_var !== nothing
        value = get(ENV, env_var, "")
        !isempty(value) && return value
    end
    path_value = Sys.which(command_name)
    path_value !== nothing && return path_value
    return fallback === nothing ? String(command_name) : String(fallback)
end

function resolve_existing_path(
    cli_value::Union{Nothing,AbstractString},
    env_var::AbstractString,
    default_path::AbstractString,
    label::AbstractString,
)
    path =
        cli_value === nothing ? get(ENV, env_var, String(default_path)) : String(cli_value)
    isfile(path) || throw(ArgumentError("$label not found: $path"))
    return path
end

struct ProcessResult
    args::Vector{String}
    code::Int
    stdout::String
    stderr::String
end

function run_checked(
    args::AbstractVector{<:AbstractString}; cwd::Union{Nothing,AbstractString}=nothing
)
    argv = String.(args)
    isempty(argv) && throw(ArgumentError("external command cannot be empty"))
    command = Cmd(argv)
    cwd === nothing || (command = Cmd(command; dir=String(cwd)))
    stdout_buffer = IOBuffer()
    stderr_buffer = IOBuffer()
    process = run(pipeline(command; stdout=stdout_buffer, stderr=stderr_buffer); wait=false)
    wait(process)
    result = ProcessResult(
        argv, process.exitcode, String(take!(stdout_buffer)), String(take!(stderr_buffer))
    )
    result.code == 0 || throw(
        ErrorException(
            "Command failed (exit code $(result.code)): $(join(result.args, " "))\n" *
            "stdout:\n$(result.stdout)\nstderr:\n$(result.stderr)",
        ),
    )
    return result
end
