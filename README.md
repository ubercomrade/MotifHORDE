# MotifHORDE

MotifHORDE is a Julia pipeline for reproducible de novo motif discovery. It
trains each requested parameter combination on odd/even foreground splits,
evaluates held-out motifs, matches reproducible models, reruns discovery on the
full foreground set, and writes ranked final models.

The runtime is `MotifHORDE.jl` plus `Mimosa.jl`. Python, pickle/joblib, C++ and
`mimosa-tool` are not runtime dependencies.

## Install

Requirements:

- Julia 1.12 or a later Julia 1.x release supported by the pinned environment;
- MEME Suite for `streme` and `meme`;
- BaMMmotif for `bamm`;
- Java and externally supplied Jstacs JARs for `dimont` and `slim`;
- the independent Julia SiteGA project for `sitega`.

Instantiate the package:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The installed `Mimosa.jl` revision is recorded in `Manifest.toml`. Its public
API is the only model reader, scanner, site/PFM extractor, serializer and
comparator used by this package.

## Usage

```bash
julia --project=. bin/motifhorde.jl \
  foreground.fa background.fa promoters.fa output/ \
  --tool streme --length 8-20-4 --nmotifs 5 --metric pauROC
```

Supported tools are `streme`, `meme`, `bamm`, `dimont`, `slim`, and `sitega`.
Ranges accept `start-end-step`, comma-separated values, or a single value.
Use `--jobs` for the shared bootstrap budget and `--seed` for deterministic
per-task seeds.

Executable resolution is explicit: CLI override, environment variable, PATH,
then the documented fallback. The variables are `HORDEMOTIFS_STREME_COMMAND`,
`HORDEMOTIFS_MEME_COMMAND`, `HORDEMOTIFS_BAMM_COMMAND`, and
`HORDEMOTIFS_SITEGA_COMMAND`. JARs use `HORDEMOTIFS_DIMONT_JAR` and
`HORDEMOTIFS_SLIM_JAR`.

## Results

Results are written below `output/<tool>/`:

```text
bootstrap/
  models/<rank>_<name>/       # Mimosa portable bundle
  statistics.json
motifs/
  models/<rank>_<name>/       # Mimosa portable bundle
  models/all_motifs_in_pfm_form.meme
  statistics.json
```

The portable model bundle is Mimosa's versioned `manifest.toml` plus binary
arrays. Old Python pickle/joblib outputs are intentionally not read or written;
regenerate them from discovery input or MEME output. See
`docs/compatibility.md` for the contract and comparison policy.

## SiteGA Contract

SiteGA is an independent process. MotifHORDE invokes the executable without
shell interpolation:

```text
sitega --foreground FG --background BG --output OUT --manifest OUT/sitega.manifest.json \
  --length LENGTH --lpd LPD --motifs COUNT [--seed SEED] [--threads JOBS]
```

The manifest must be JSON with `schema_version: 1`, `status: "success"`, and
an ordered `models` array. Every model entry contains `model_file`,
`model_type`, and `length`, with an optional stable `name`. Model files must be
readable by Mimosa. A missing, invalid, failed, or incompatible manifest is a
hard error.

## Development

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The package tests use fake external processes for the SiteGA contract. Real
data smoke tests remain opt-in and require the external tools to be installed.
