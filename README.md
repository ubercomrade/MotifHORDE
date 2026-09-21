# MotifHORDE

MotifHORDE orchestrates de novo motif discovery, odd/even bootstrap
validation, Mimosa profile comparison, and portable result export.

## Requirements

- Python `>=3.12,<3.14`;
- `mimosa-tool>=1.4.6,<1.5`;
- external discovery programs on `PATH` or supplied explicitly.

The conda-first environment includes MEME Suite and `bammmotif2`:

```bash
conda env create -f environment.yml
conda activate motifhorde
uv sync --locked
```

Dimont, SlimDimont, and SiteGA are external processes. MotifHORDE does not
ship their JARs, native code, or model artifacts.

## CLI

The four positional arguments are stable:

```text
motifhorde FOREGROUND.fa BACKGROUND.fa PROMOTERS.fa OUTPUT_DIR
```

Supported tools are `streme`, `meme`, `bamm`, `dimont`, `slim`, and `sitega`.
For example:

```bash
motifhorde peaks.fa background.fa promoters.fa output/ \
  --tool streme --length 8-20-4 --nmotifs 5 \
  --mimosa-metric cosine --jobs 4 --seed 42
```

Only Mimosa profile metrics `co`, `dice`, and `cosine` are supported. Use
`--comparison-criterion p-value` with a Mimosa null-distribution directory
bundle passed to `--mimosa-null-distribution`.

SiteGA uses `--lpd` for the maximum LPD distance (`1..6`) and `--features` for
the feature-count grid:

```bash
motifhorde peaks.fa background.fa promoters.fa output/ \
  --tool sitega --length 10-16-2 --lpd 6 --features 10-40-10 \
  --sitega-command sitega
```

Executable resolution is CLI option, environment variable, then `PATH`:

| Tool | CLI option | Environment variable |
| --- | --- | --- |
| MEME | `--meme-command` | `HORDEMOTIFS_MEME_COMMAND` |
| STREME | `--streme-command` | `HORDEMOTIFS_STREME_COMMAND` |
| BaMM | `--bamm-command` | `HORDEMOTIFS_BAMM_COMMAND` |
| SiteGA | `--sitega-command` | `HORDEMOTIFS_SITEGA_COMMAND` |
| Dimont JAR | `--dimont-jar` | `HORDEMOTIFS_DIMONT_JAR` |
| SlimDimont JAR | `--slim-jar` | `HORDEMOTIFS_SLIM_JAR` |

## Output

Results are written below `OUTPUT_DIR/<tool>`:

```text
bootstrap/
  tasks/task_0000/...
  models/<task-id>
  statistics.json
motifs/
  models/rank_001_<name>
  models/all_motifs_in_pfm_form.meme
  statistics.json
```

Every model is a Mimosa portable directory bundle. `statistics.json` contains
`schema_version`, `tool`, and an ordered `motifs` list. Bootstrap entries have
`id`, `model`, `name`, `split`, `params`, and `metrics`; final entries have
`rank`, `model`, `name`, `params`, and `metrics`.

Existing result bundles are never overwritten. Regenerate old pickle/joblib
outputs from discovery inputs or interoperable MEME data.

## SiteGA manifest

The external SiteGA process receives `--manifest <path>` and writes a versioned
JSON manifest. MotifHORDE validates its schema/status, ordered model list,
relative paths, containment, file existence, `sitega` type, and model length.
It does not infer model order from filenames.

## Python API

Use Mimosa directly for model and sequence operations:

```python
from mimosa import BestPerSequence, TopFractionHits, read_fasta, read_model
from mimosa import reconstruct_pfm, scan

sequences, names = read_fasta("peaks.fa")
model = read_model("motif.meme", format="meme")
scores = scan(model, sequences, strands="best")
pfm = reconstruct_pfm(model, sequences, TopFractionHits(0.10, BestPerSequence()))
```

MotifHORDE-specific helpers are limited to external discovery adapters,
evaluation, orchestration, discovery FASTA variants, and combined MEME PFM
export. There is no `GenericModel`, pickle loader, Tomtom comparator, or
embedded SiteGA implementation.

## Development checks

```bash
uv lock
uv sync --locked
uv run pytest
uv run ruff check .
uv run ruff format --check .
uv build
git diff --check
```
