# MotifHORDE

**MotifHORDE** stands for **Motif Harmonized Orchestration for Robust Discovery
and Evaluation**.

`motifhorde` runs de novo motif discovery with odd/even bootstrap validation,
motif comparison, evaluation, and final model selection. It harmonizes outputs
from multiple motif discovery tools into one consistent Python API and pipeline.

The runtime model API uses the `mimosa` style:

- `GenericModel` stores `type_key`, `name`, `representation`, `length`, and `config`.
- Model operations are functional: `read_model`, `scan_model`, `get_pfm`, `get_sites`, `calculate_threshold_table`.
- Sequence inputs are dense batches: `{"values", "lengths", "padding_value"}`.
- Score/profile outputs are dense masked batches or profile bundles.

Legacy model classes and ragged payloads are removed from production code.

## Installation

For command-line use with MEME/STREME, create a conda or mamba environment:

```bash
mamba env create -f environment.yml
conda activate motifhorde
```

The environment installs MEME Suite from Bioconda, so `meme` and `streme` are
available through `PATH`. Dimont and SlimDimont JAR files are bundled with the
Python package.

## Python API

```python
from motifhorde.io import read_fasta
from motifhorde.models import read_model, scan_model, get_pfm, get_sites

sequences = read_fasta("peaks.fa")
model = read_model("motif.meme", "pwm")

scores = scan_model(model, sequences, strand="best")
pfm = get_pfm(model, sequences, top_fraction=0.10)
sites = get_sites(model, sequences)
```

## CLI

```bash
motifhorde peaks.fa background.fa promoters.fa output/ -t streme -l 8-20-4
motifhorde peaks.fa background.fa promoters.fa output/ -t bamm -l 10-14-2 -o 1-4-1
motifhorde peaks.fa background.fa promoters.fa output/ -t sitega -l 10-16-2 --lpd 10-40-10
```

Supported discovery tools are `streme`, `meme`, `bamm`, `dimont`, `slim`, and `sitega`.
Supported comparators are `tomtom` and `continuous`.
Continuous profile metrics are `co`, `co_rowwise`, `dice`, `dice_rowwise`, and `cosine`.

BaMM models currently use the uniform-background reader inherited from the
`mimosa` model implementation.
