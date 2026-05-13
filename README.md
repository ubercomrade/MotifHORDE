# MotifHORDE

**MotifHORDE** stands for **Motif Harmonized Orchestration for Robust Discovery
and Evaluation**.

`motifhorde` is a command-line and Python toolkit for de novo transcription
factor motif discovery, validation, comparison, and final model selection.

The core goal is practical parameter selection for motif discovery tools. For a
grid of discovery parameters, MotifHORDE trains motifs on odd/even splits of
ranked ChIP-seq peaks, evaluates them on the held-out split, matches similar
motifs between validation folds, and then reruns discovery on the full dataset
with selected parameters.

## Method Overview

MotifHORDE implements a parameter-selection workflow for de novo motif discovery
on ranked TF-bound regions. The input foreground FASTA is expected to preserve
the peak order from the upstream ChIP-seq peak caller or ranking procedure. This
matches common ChIP-seq workflows where peaks are ordered by signal strength,
for example MACS peak calling and ranking
([Zhang et al., 2008](https://doi.org/10.1186/gb-2008-9-9-r137)).

The method is designed for comparing motif models learned with different
discovery parameters. For PWM tools, the primary grid parameter is usually motif
length. For Markov and dependency-aware models, the grid can also include model
order or locally positioned dinucleotide settings. This reflects the broader
motif-discovery setting where PWM models remain a common baseline
([Machanick & Bailey, 2011](https://doi.org/10.1093/bioinformatics/btr189);
[Bailey, 2021](https://doi.org/10.1093/bioinformatics/btab203)), while
alternative models can capture dependencies that PWMs miss
([Grau et al., 2013](https://doi.org/10.1093/nar/gkt831);
[Siebert & Soding, 2016](https://doi.org/10.1093/nar/gkw521);
[Tsukanov et al., 2021](https://doi.org/10.18699/VJ21.002);
[Tsukanov et al., 2022](https://doi.org/10.3389/FPLS.2022.938545)).

The implemented workflow is:

1. **Build the parameter grid.** The CLI expands ranges such as `8-20-4`,
   comma-separated values such as `8,10,12`, or a single value. The pipeline
   then evaluates the Cartesian product of all tool-specific parameters.
2. **Split foreground peaks by odd/even index.** The split uses 1-based peak
   indices and preserves the original signal-ranked order.
3. **Run two-fold validation for each parameter combination.** On the first
   fold, motifs are learned from odd peaks and evaluated on even peaks. On the
   second fold, motifs are learned from even peaks and evaluated on odd peaks.
4. **Evaluate each validation motif.** Positives are the held-out foreground
   peaks. Negatives are the provided background FASTA. Motif scores are converted
   into binary-classification curves and summarized as AUROC, AUPRC, partial
   AUROC, and partial AUPRC.
5. **Select nonredundant validation motifs within each parameter set.** Odd-fold
   and even-fold motifs are reduced separately, using the selected validation
   metric to keep the best representative among similar motifs.
6. **Compare selected odd-fold and even-fold motifs with identical parameters.**
   This step identifies motifs that are reproducible across the two validation
   folds without mixing different parameter settings. MotifHORDE supports a
   Tomtom-like matrix comparator, following the
   alignment-based idea of Tomtom
   ([Gupta et al., 2007](https://doi.org/10.1186/gb-2007-8-2-r24)), and a
   profile comparator that compares recognition profiles, following the idea
   that models can be compared by predicted binding or affinity patterns rather
   than only by matrix columns
   ([Vorontsov et al., 2013](https://doi.org/10.1186/1748-7188-8-23);
   [Lambert et al., 2016](https://doi.org/10.1093/bioinformatics/btw489)).
7. **Deduplicate matched validation pairs within each parameter set.**
   Comparison results are sorted by p-value when Monte Carlo p-values are
   available, otherwise by similarity score. The pipeline then keeps one match
   per query motif and one match per target motif for the current parameter set.
8. **Assign validation performance to matched pairs.** For every retained
   odd/even pair, MotifHORDE averages the validation metrics of the two matched
   motifs. These averaged metrics represent how well a reproducible motif
   performed under the tested parameter combination.
9. **Rerun discovery on the full foreground set.** For every retained parameter
   group, the selected discovery tool is run again on all foreground peaks.
10. **Select final full-data motifs.** Each full-data motif is compared with the
    corresponding odd/even validation references. A full-data motif is retained
    only if it is similar to both references, then the best match is selected
    using the comparator direction: lower averaged p-value for p-value based
    comparisons, or higher averaged score for score based comparisons.
11. **Globally remove redundant final motifs.** After all retained parameter
    groups have been processed, selected full-data motifs are sorted by the
    requested validation metric and redundant motifs across parameter settings
    are removed.
12. **Export final results.** Final motifs are sorted by the requested selection
    metric (`--metric`, default `pauROC`). Pickled `GenericModel` objects are
    written to the output directory, and all final motifs are also exported as a
    MEME-format PFM file. For non-PWM models, the PFM is reconstructed from
    predicted sites on the provided promoter sequences.

Across the selection pipeline, motifs are treated as similar when comparison
results have `p-value <= 0.001` or `score >= 0.9`. These method-level thresholds
are fixed inside the pipeline rather than exposed as separate selection flags.

## Supported Models

MotifHORDE uses a functional runtime API inspired by `mimosa`.

- `GenericModel` stores `type_key`, `name`, `representation`, `length`, and
  `config`.
- Model operations are functions: `read_model`, `write_model`, `scan_model`,
  `scan_model_strands`, `get_pfm`, `get_sites`, and
  `calculate_threshold_table`.
- Sequence inputs are dense batches with `values`, `lengths`, and
  `padding_value`.
- Score/profile outputs are dense masked batches or strand-aware profile
  bundles.

Supported model readers and scanners:

| Type key | Main file formats | Notes |
| :--- | :--- | :--- |
| `pwm` | `.meme`, `.txt`, `.pfm`, `.pkl` | MEME/PFM matrices are converted to PWM log-odds. |
| `bamm` | `.ihbcp` | Read as a dense Markov tensor with uniform-background log-odds. |
| `sitega` | `.mat`, `.pkl` | Locally positioned dinucleotide model. |
| `dimont` | `.xml`, `.pkl` | Jstacs Dimont XML is converted to a dense scoring tensor. |
| `slim` | `.xml`, `.pkl` | Jstacs SlimDimont XML is converted to a dense scoring tensor. |
| `scores` | FASTA-like numeric profiles | Used for precomputed score/profile workflows. |

Legacy motif classes and ragged payloads are not used in production code.

## Discovery Tools

The CLI can run the following external discovery tools:

| CLI value | Model family | Main parameter grid |
| :--- | :--- | :--- |
| `streme` | PWM | `--length` |
| `meme` | PWM | `--length` |
| `bamm` | BaMM | `--length`, `--order` |
| `dimont` | Dimont | `--length` plus Jstacs options |
| `slim` | SlimDimont | `--length` plus Jstacs options |
| `sitega` | SiteGA | `--length`, `--lpd` |

Ranges use one of three formats:

- `start-end-step`, for example `8-20-2`;
- comma-separated values, for example `8,10,12`;
- a single value, for example `12`.

## Motif Comparison

MotifHORDE includes two pipeline comparators.

### Matrix Comparator

`-c tomtom` uses a Tomtom-like alignment over motif matrices or tensors. It
supports:

- `--tomtom-metric pcc` for Pearson correlation;
- `--tomtom-metric ed` for negative Euclidean distance;
- column permutation-based Monte Carlo null scores with `--tomtom-perm`;
- optional row permutation with `--tomtom-permute-rows`;
- optional PFM reconstruction with `--pfm-mode`.

When `--pfm-mode` is enabled, or when two models have different type keys, the
comparison derives PFMs by scanning sequences and using the top-scoring predicted
sites. This makes heterogeneous model comparison possible without requiring both
models to share the same internal representation.

The CLI exposes `--tomtom-pval` for compatibility with existing invocations, but
pipeline selection uses the fixed method-level similarity rule: `p-value <=
0.001` when p-values are available, otherwise `score >= 0.9`. Tomtom-like
results are sorted by p-value when Monte Carlo p-values are present.

### Continuous Profile Comparator

`-c continuous` compares the recognition profiles produced by motif models on
the same sequence set. Models are compared by their functional output rather
than only by their internal parameters.

Supported metrics:

- `co`: continuous overlap coefficient;
- `co_rowwise`: row-wise continuous overlap;
- `dice`: continuous Dice coefficient;
- `dice_rowwise`: row-wise continuous Dice;
- `cosine`: row-wise cosine similarity.

Profile comparison supports Monte Carlo null estimation through surrogate
profile generation. The surrogate generator applies randomized convolutional
distortion controlled by `--c-distortion`, `--min-kernel-size`, and
`--max-kernel-size`.

## Evaluation

Each bootstrap motif is evaluated as a binary classifier:

- positives are held-out foreground peaks;
- negatives are the provided background FASTA;
- foreground scores use the best site per sequence;
- background scores use either all sites or the best site per sequence,
  controlled by `--background-type sites|peaks`.

Reported metrics:

- `auROC`;
- `auPRC`;
- `pauROC`;
- `pauPRC`.

Partial metrics are calculated at the configured false-positive-rate threshold
(`--fpr`, default `0.001`). Final motifs are ranked by `--metric`, which defaults
to `pauROC`.

## Installation

MotifHORDE requires Python 3.12 or newer.

For command-line use with MEME/STREME, create the provided conda or mamba
environment:

```bash
conda env create -f environment.yml
conda activate motifhorde
```

The environment installs MEME Suite from Bioconda, so `meme` and `streme` are
available through `PATH`.

For local development from this repository:

```bash
uv sync
uv pip install -e .
```

External dependency notes:

- `streme` and `meme` require MEME Suite.
- `bamm` requires `BaMMmotif` and uses STREME for PWM initialization.
- `dimont` and `slim` require Java; `Dimont.jar` and `SlimDimont.jar` are
  bundled in the Python package.
- `sitega` requires `andy05cell.exe` to be available through `PATH`.

Executable and JAR paths can be overridden by CLI flags or environment variables:

| Dependency | CLI flag | Environment variable |
| :--- | :--- | :--- |
| MEME | `--meme-command` | `HORDEMOTIFS_MEME_COMMAND` |
| STREME | `--streme-command` | `HORDEMOTIFS_STREME_COMMAND` |
| BaMMmotif | `--bamm-command` | `HORDEMOTIFS_BAMM_COMMAND` |
| Dimont JAR | `--dimont-jar` | `HORDEMOTIFS_DIMONT_JAR` |
| SlimDimont JAR | `--slim-jar` | `HORDEMOTIFS_SLIM_JAR` |

## CLI Usage

```bash
motifhorde foreground.fa background.fa promoters.fa output/ \
  -t streme \
  -l 8-20-4 \
  -n 5 \
  -m pauROC
```

The required positional inputs are:

| Argument | Meaning |
| :--- | :--- |
| `foreground` | FASTA with ranked foreground peak sequences. |
| `background` | FASTA with negative/background sequences. |
| `promoters` | FASTA used to derive final PFM representations for selected motifs. |
| `output` | Output directory. |

Examples:

```bash
# PWM discovery with STREME over motif lengths 8, 12, 16, 20
motifhorde peaks.fa background.fa promoters.fa output/ \
  -t streme \
  -l 8-20-4

# BaMM discovery over motif length and Markov order grids
motifhorde peaks.fa background.fa promoters.fa output/ \
  -t bamm \
  -l 10-14-2 \
  -o 1-4-1

# SiteGA discovery over motif length and LPD count grids
motifhorde peaks.fa background.fa promoters.fa output/ \
  -t sitega \
  -l 10-16-2 \
  --lpd 10-40-10

# Continuous profile comparison instead of matrix comparison
motifhorde peaks.fa background.fa promoters.fa output/ \
  -t dimont \
  -l 10,12,14 \
  -c continuous \
  --c-metric co \
  --c-perm 1000

# Jstacs tools with explicit Java settings
motifhorde peaks.fa background.fa promoters.fa output/ \
  -t slim \
  -l 12 \
  --java-xmx 8G \
  --jstacs-threads 4
```

## Output Layout

For a selected tool, results are written under:

```text
output/
  <tool>/
    bootstrap/
      statistics.json
      models/
        0000_<motif>.pkl
        ...
    motifs/
      statistics.json
      models/
        all_motifs_in_pfm_form.meme
        001_<motif>.pkl
        ...
```

`bootstrap/statistics.json` stores metrics for odd/even validation motifs.
`motifs/statistics.json` stores metrics associated with selected final motifs.
Final models are saved as pickled `GenericModel` objects. A MEME-format file
with all selected final motifs in PFM form is also written for compatibility
with downstream PWM-oriented tools.

## Python API

```python
from motifhorde.io import read_fasta
from motifhorde.models import (
    calculate_threshold_table,
    get_pfm,
    get_sites,
    read_model,
    scan_model,
)

sequences = read_fasta("peaks.fa")
background = read_fasta("background.fa")
model = read_model("motif.meme", "pwm")

scores = scan_model(model, sequences, strand="best")
threshold_table = calculate_threshold_table(model, background, strand="best")
pfm = get_pfm(model, sequences, top_fraction=0.10)
sites = get_sites(
    model,
    sequences,
    mode="threshold",
    fpr_threshold=0.001,
    background_sequences=background,
    threshold_table=threshold_table,
)
```

Direct comparison API:

```python
from motifhorde.comparison import TomtomComparator, UniversalMotifComparator

tomtom = TomtomComparator(metric="pcc", n_permutations=1000, seed=1)
matrix_results = tomtom.compare([model_a], [model_b], sequences=sequences)

continuous = UniversalMotifComparator(
    metric="co",
    n_permutations=1000,
    filter_type="p-value",
    filter_threshold=0.05,
    seed=1,
)
profile_results = continuous.compare([model_a], [model_b], sequences=sequences)
```

Programmatic full pipeline:

```python
from motifhorde.comparison import TomtomComparator
from motifhorde.discovery import StremeDiscoveryTool
from motifhorde.evaluation import PerformanceEvaluator
from motifhorde.pipeline import DeNovoPipeline

pipeline = DeNovoPipeline(
    discovery_tool=StremeDiscoveryTool(),
    evaluator=PerformanceEvaluator(background_type="peaks"),
    comparator=TomtomComparator(metric="pcc", n_permutations=1000),
    fpr_threshold=0.001,
    number_of_motifs=5,
)

pipeline.run(
    "peaks.fa",
    "background.fa",
    "promoters.fa",
    "output/",
    discovery_params={"length": [8, 12, 16, 20]},
    metric="pauROC",
)
```

## Bibliography

Bailey, T. L. (2021). STREME: accurate and versatile sequence motif discovery.
Bioinformatics, 37(18), 2834-2840. https://doi.org/10.1093/bioinformatics/btab203

Grau, J., Posch, S., Grosse, I., & Keilwagen, J. (2013). A general approach for
discriminative de novo motif discovery from high-throughput data. Nucleic Acids
Research, 41(21), e197. https://doi.org/10.1093/nar/gkt831

Gupta, S., Stamatoyannopoulos, J. A., Bailey, T. L., & Noble, W. (2007).
Quantifying similarity between motifs. Genome Biology, 8(2), R24.
https://doi.org/10.1186/gb-2007-8-2-r24

Kulakovskiy, I. V., Levitsky, V., Oshchepkov, D., Bryzgalov, L., Vorontsov, I.,
& Makeev, V. (2013). From binding motifs in ChIP-Seq data to improved models of
transcription factor binding sites. Journal of Bioinformatics and Computational
Biology, 11(1), 1340004. https://doi.org/10.1142/S0219720013400040

Lambert, S. A., Albu, M., Hughes, T. R., & Najafabadi, H. S. (2016). Motif
comparison based on similarity of binding affinity profiles. Bioinformatics,
32(22), 3504-3506. https://doi.org/10.1093/bioinformatics/btw489

Machanick, P., & Bailey, T. L. (2011). MEME-ChIP: motif analysis of large DNA
datasets. Bioinformatics, 27(12), 1696-1697.
https://doi.org/10.1093/bioinformatics/btr189

Mahony, S., & Benos, P. V. (2007). STAMP: A web tool for exploring DNA-binding
motif similarities. Nucleic Acids Research, 35(Web Server issue), W253-W258.
https://doi.org/10.1093/nar/gkm272

Siebert, M., & Soding, J. (2016). Bayesian Markov models consistently outperform
PWMs at predicting motifs in nucleotide sequences. Nucleic Acids Research,
44(13), 6055-6069. https://doi.org/10.1093/nar/gkw521

Tsukanov, A. V., Levitsky, V. G., & Merkulova, T. I. (2021). Application of
alternative de novo motif recognition models for analysis of structural
heterogeneity of transcription factor binding sites: A case study of FOXA2
binding sites. Vavilov Journal of Genetics and Breeding, 25(1), 7-17.
https://doi.org/10.18699/VJ21.002

Tsukanov, A. V., Mironova, V. V., & Levitsky, V. G. (2022). Motif models
proposing independent and interdependent impacts of nucleotides are related to
high and low affinity transcription factor binding sites in Arabidopsis.
Frontiers in Plant Science, 13, 2637.
https://doi.org/10.3389/FPLS.2022.938545

Vorontsov, I. E., Kulakovskiy, I. V., & Makeev, V. J. (2013). Jaccard index
based similarity measure to compare transcription factor binding site models.
Algorithms for Molecular Biology, 8(1), 23.
https://doi.org/10.1186/1748-7188-8-23

Zhang, Y., Liu, T., Meyer, C. A., Eeckhoute, J., Johnson, D. S., Bernstein,
B. E., Nusbaum, C., Myers, R. M., Brown, M., Li, W., & Liu, X. S. (2008).
Model-based analysis of ChIP-Seq (MACS). Genome Biology, 9(9), R137.
https://doi.org/10.1186/gb-2008-9-9-r137
