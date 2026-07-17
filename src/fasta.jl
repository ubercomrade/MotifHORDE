"""FASTA and interoperable model output at the I/O boundary."""

function read_fasta(path::AbstractString; max_sequences::Integer=1_000_000)
    isfile(path) || throw(ArgumentError("FASTA file not found: $path"))
    filesize(path) == 0 && return (Mimosa.empty_sequence_batch(), String[])
    return Mimosa.read_fasta(path; max_sequences=Int(max_sequences))
end

function _base_string(sequence::AbstractVector{UInt8})
    decoder = UInt8['A', 'C', 'G', 'T', 'N']
    return String(UInt8[decoder[min(Int(code) + 1, 5)] for code in sequence])
end

function write_fasta(batch::Mimosa.EncodedSequenceBatch, path::AbstractString)
    open(path, "w") do io
        for index in 1:Mimosa.nsequences(batch)
            println(io, ">seq$(index - 1)")
            println(io, _base_string(Mimosa.sequence(batch, index)))
        end
    end
    return path
end

function write_jstacs_fasta(
    input_path::AbstractString,
    output_path::AbstractString;
    position_tag::AbstractString="position",
    value_tag::AbstractString="value",
)
    if !isfile(input_path) || filesize(input_path) == 0
        open(output_path, "w") do _
        end
        return output_path
    end
    batch, _ = read_fasta(input_path)
    open(output_path, "w") do io
        for index in 1:Mimosa.nsequences(batch)
            length = Mimosa.seqlength(batch, index)
            println(
                io,
                "> $(position_tag): $(length ÷ 2); $(value_tag): $(@sprintf("%.1f", index))",
            )
            println(io, _base_string(Mimosa.sequence(batch, index)))
        end
    end
    return output_path
end

function write_meme(
    pfms::AbstractVector, info::AbstractVector{<:Tuple}, path::AbstractString
)
    length(pfms) == length(info) || throw(ArgumentError("PFM metadata length mismatch."))
    open(path, "w") do io
        println(io, "MEME version 4\n")
        println(io, "ALPHABET= ACGT\n")
        println(io, "strands: + -\n")
        println(io, "Background letter frequencies")
        println(io, "A 0.25 C 0.25 G 0.25 T 0.25\n")
        for (pfm, metadata) in zip(pfms, info)
            name, width = metadata
            size(pfm, 1) == 4 || throw(ArgumentError("PFM for $name must have 4 rows."))
            println(io, "MOTIF $name")
            println(io, "letter-probability matrix: alength= 4 w= $width nsites= 20 E= 0")
            for position in 1:size(pfm, 2)
                println(
                    io,
                    join((@sprintf(" %.6f", Float64(pfm[base, position])) for base in 1:4)),
                )
            end
            println(io)
        end
    end
    return path
end
