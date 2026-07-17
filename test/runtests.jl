using Test
using Aqua
using ArgParse
using Mimosa
using MotifHORDE

@testset "MotifHORDE" begin
    @test MotifHORDE.parse_range("8-20-4") == [8, 12, 16, 20]
    @test MotifHORDE.parse_range("8, 10, 12") == [8, 10, 12]
    @test_throws ArgumentError MotifHORDE.parse_range("8-20-0")

    @testset "CLI parsing" begin
        parsed = ArgParse.parse_args(
            ["foreground.fa", "background.fa", "promoters.fa", "output"],
            MotifHORDE._cli_settings(),
        )
        @test parsed["tool"] == "streme"
        @test parsed["nmotifs"] == 5
        @test parsed["lpd"] == 6
        @test parsed["features"] == "10-40-5"
        @test parsed["verbose"] == false
        @test_throws ArgumentError MotifHORDE.parse_range("8-20-0")
    end

    @testset "FASTA and Jstacs I/O" begin
        mktempdir() do directory
            input = joinpath(directory, "input.fa")
            annotated = joinpath(directory, "annotated.fa")
            empty_input = joinpath(directory, "empty.fa")
            write(empty_input, "")
            empty_batch, empty_names = MotifHORDE.read_fasta(empty_input)
            @test Mimosa.nsequences(empty_batch) == 0
            @test isempty(empty_names)
            write(input, ">a\nACGTAC\n>b\nTTTTAA\n")
            MotifHORDE.write_jstacs_fasta(input, annotated)
            @test read(annotated, String) ==
                "> position: 3; value: 1.0\nACGTAC\n> position: 3; value: 2.0\nTTTTAA\n"
            batch, names = MotifHORDE.read_fasta(input)
            @test names == ["a", "b"]
            output = joinpath(directory, "roundtrip.fa")
            MotifHORDE.write_fasta(batch, output)
            reread, _ = MotifHORDE.read_fasta(output)
            @test [Mimosa.sequence(batch, i) for i in 1:Mimosa.nsequences(batch)] == [Mimosa.sequence(reread, i) for i in 1:Mimosa.nsequences(reread)]
        end
    end

    @testset "Model API and evaluation" begin
        pfm = Float32[
            0.8 0.2 0.2 0.2;
            0.1 0.3 0.3 0.3;
            0.05 0.25 0.25 0.25;
            0.05 0.25 0.25 0.25
        ]
        model = Mimosa.pwm_from_pfm(pfm; name="m1")
        sequences = Mimosa.EncodedSequenceBatch([
            Mimosa.encode_sequence("ACGTACGT"), Mimosa.encode_sequence("TTTTTTTT")
        ])
        @test length(MotifHORDE.best_scores(model, sequences)) == 2
        @test haskey(
            MotifHORDE.evaluate(
                MotifHORDE.PerformanceEvaluator(background_type=:peaks),
                model,
                sequences,
                sequences,
                0.1,
            ),
            "auROC",
        )
        mktempdir() do directory
            path = joinpath(directory, "model")
            MotifHORDE.write_model(path, model)
            loaded = MotifHORDE.read_model(path)
            @test loaded == model
        end
    end

    @testset "Comparison rules" begin
        @test MotifHORDE.default_threshold_for_criterion("score") == 0.9
        @test MotifHORDE.default_threshold_for_criterion("p-value") == 0.05
        @test MotifHORDE.comparison_column_for_criterion("p-value") == Symbol("adj.p-value")
        @test MotifHORDE.format_params(Dict("z" => 1, "a" => 2)) == "a-2,z-1"
    end

    @testset "SiteGA process contract" begin
        mktempdir() do directory
            foreground = joinpath(directory, "foreground.fa")
            background = joinpath(directory, "background.fa")
            output = joinpath(directory, "output")
            fake = joinpath(directory, "sitega")
            write(foreground, ">x\nACGTACGTACGT\n>y\nTTTTACGTACGT\n")
            write(background, ">x\nTTTTCCCCGGGG\n>y\nGGGGAAAACCCC\n")
            mkpath(output)
            model = Mimosa.pwm_from_pfm(
                Float32[
                    0.8 0.2 0.2 0.2;
                    0.1 0.3 0.3 0.3;
                    0.05 0.25 0.25 0.25;
                    0.05 0.25 0.25 0.25
                ];
                name="external",
            )
            MotifHORDE.write_model(joinpath(output, "model"), model)
            write(
                fake,
                raw"""#!/bin/sh
                set -eu
                output=""
                manifest=""
                while [ "$#" -gt 0 ]; do
                    case "$1" in
                        --output) output="$2"; shift 2 ;;
                        --manifest) manifest="$2"; shift 2 ;;
                        *) shift ;;
                    esac
                done
                printf '%s\n' '{"schema_version":1,"status":"success","message":"ok","models":[{"model_file":"model","model_type":"pwm","length":4,"name":"Sitega-1"}]}' > "$manifest"
                """,
            )
            chmod(fake, 0o755)
            models = MotifHORDE.discover(
                MotifHORDE.SitegaDiscoveryTool(command=fake, threads=1, seed=7),
                foreground,
                background,
                output,
                1;
                length=4,
                lpd=3,
            )
            @test [
                (MotifHORDE.motif_name(model), MotifHORDE.motif_length(model)) for
                model in models
            ] == [("Sitega-1", 4)]
        end
    end
end

@testset "Aqua" begin
    Aqua.test_all(MotifHORDE; deps_compat=false, stale_deps=false)
end
