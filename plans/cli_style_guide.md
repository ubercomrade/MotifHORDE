# CLI Standardization Plan for `hordemotifs/cli.py`

## 1. Style Guide

### Capitalization & Punctuation
- **Help Strings:** Start with a capital letter. Do not end with a period unless the help string consists of multiple sentences.
- **Argument Names:** Use lowercase with hyphens (kebab-case) for long options (e.g., `--background-type`).
- **Defaults:** Append `(default: %(default)s)` to the end of help strings where a default value is provided.

### Terminology
- **Files:** Use specific terms like "foreground FASTA file" instead of "foreground sequences".
- **Processes:** Use "significance testing" instead of just "p-value calculation".
- **Parallelism:** Use "parallel jobs" consistently.

## 2. Revised Text Strings

### Command Description
**Old:** `HordeMotifs: De novo motif discovery with odd/even bootstrap validation`
**New:** `HordeMotifs: De novo motif discovery pipeline with odd/even bootstrap validation`

### Argument Help Strings Mapping

| Argument | New Help String |
| :--- | :--- |
| `foreground` | Path to the foreground FASTA file containing peak sequences |
| `background` | Path to the background FASTA file |
| `promoters` | Path to the promoter FASTA file used for threshold calculation |
| `output` | Path to the output directory where results will be saved |
| `--tool` | De novo motif discovery tool to use (default: %(default)s) |
| `--nmotifs` | Number of motifs to discover per run (default: %(default)s) |
| `--length` | Range of motif lengths to discover. Format: 'start-end-step' (e.g., 8-20-4), comma-separated list (e.g., 8,10,12), or a single value (default: %(default)s) |
| `--order` | Range of Markov model orders. Format: 'start-end-step' (e.g., 1-4-1), comma-separated list (e.g., 1,2,3), or a single value (default: %(default)s) |
| `--lpd` | Range of locally positioned dinucleotide (LPD) distances for SiteGA. Format: 'start-end-step', comma-separated list, or a single value (default: %(default)s) |
| `--fpr` | False Positive Rate (FPR) threshold for partial AUC calculation (default: %(default)s) |
| `--background-type` | Method for background scoring. 'peaks' uses the best site per sequence; 'sites' uses all sites (default: %(default)s) |
| `--metric` | Performance metric used to select the best motifs (default: %(default)s) |
| `--comparator` | Method used for comparing discovered motifs (default: %(default)s) |
| `--tomtom-metric` | Distance metric for TomTom comparison. 'pcc' is Pearson Correlation Coefficient; 'ed' is Euclidean Distance (default: %(default)s) |
| `--tomtom-pval` | P-value threshold for TomTom Monte Carlo significance (default: %(default)s) |
| `--tomtom-perm` | Number of permutations for TomTom significance testing (default: %(default)s) |
| `--tomtom-permute-rows` | Permute rows instead of columns when generating the null distribution |
| `--tomtom-jobs` | Number of parallel jobs for TomTom calculations (default: %(default)s) |
| `--cj-metric` | Metric for Jaccard/Overlap comparison. 'cj' is Jaccard index; 'co' is overlap coefficient (default: %(default)s) |
| `--cj-perm` | Number of permutations for significance testing (default: %(default)s) |
| `--cj-distortion` | Distortion level applied during surrogate motif generation (default: %(default)s) |
| `--cj-filter` | Criterion for filtering comparison results (default: %(default)s) |
| `--cj-threshold` | Numerical threshold for the filtering criterion (default: %(default)s) |
| `--cj-jobs` | Number of parallel jobs for comparisons. Set to -1 to use all available cores (default: %(default)s) |
| `--motali-threshold` | Similarity threshold for Motali comparison (default: %(default)s) |
| `--seed` | Random seed for reproducible results (default: %(default)s) |
| `--verbose` | Enable verbose output logging |

### Epilog
```python
epilog="""
Examples:
  # Basic PWM discovery with STREME
  hordemotifs peaks.fa background.fa promoters.fa output/ -t streme -l 8-20-4

  # Markov model-based motifs with different orders
  hordemotifs peaks.fa bg.fa promoters.fa output/ -t bamm -l 10-14-2 -o 1-4-1

  # SiteGA with custom LPD range
  hordemotifs peaks.fa bg.fa promoters.fa output/ -t sitega -l 10-16-2 --lpd 10-40-10

  # Single length value (when testing specific length)
  hordemotifs peaks.fa bg.fa promoters.fa output/ -t streme -l 12

  # Multiple specific values using comma-separated format
  hordemotifs peaks.fa bg.fa promoters.fa output/ -t bamm -l 10,12,14 -o 1,2,3

  # Custom comparison with Jaccard metric
  hordemotifs peaks.fa bg.fa promoters.fa output/ -c jaccard --cj-metric cj --cj-perm 5000

  # TomTom comparison with custom threshold
  hordemotifs peaks.fa bg.fa promoters.fa output/ -c tomtom --tomtom-pval 0.0001 --tomtom-metric pcc
"""
```

## 3. Standardized Error Messages

| Context | Template | Example |
| :--- | :--- | :--- |
| **File Not Found** | `ERROR: <Description> file not found: <Path>` | `ERROR: Foreground FASTA file not found: data/peaks.fa` |
| **Dependency Missing** | `ERROR: <Tool> dependency missing. Please install <Tool> from <URL>` | `ERROR: STREME dependency missing. Please install MEME Suite from https://meme-suite.org/` |
| **Invalid Argument** | `ERROR: Invalid format for <Argument>: <Value>. Expected <Format>` | `ERROR: Invalid format for --length: 8-20. Expected 'start-end-step' or comma-separated values` |
| **General Error** | `ERROR: <Message>` | `ERROR: Unknown discovery tool: magic_tool` |

## 4. Implementation Plan
1.  Update `create_arg_parser` in `hordemotifs/cli.py` with the new description, epilog, and help strings.
2.  Update `check_dependencies` to use the new error message format.
3.  Update `parse_range` and `setup_discovery_params` to use the new error message format.
4.  Update `main_cli` to use the new error message format for file checks.
