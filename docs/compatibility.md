# Compatibility Contract

## Version

Contract version: `1`.

| Area | Contract |
| --- | --- |
| Julia | Julia 1.12 environment recorded in `Manifest.toml` |
| Model runtime | Mimosa.jl `0.1.0`, pinned git revision |
| Discovery tools | `streme`, `meme`, `bamm`, `dimont`, `slim`, `sitega` |
| Output model format | Mimosa portable directory bundle, `MODEL_FORMAT_VERSION == 2` |
| Interoperable PWM output | MEME PFM file `all_motifs_in_pfm_form.meme` |
| Legacy model output | Python pickle/joblib is unsupported |

## CLI Matrix

All commands use four positional arguments: foreground FASTA, background FASTA,
promoters FASTA, and output directory.

| Option | Values/default | Validation |
| --- | --- | --- |
| `--tool` | `streme` | one supported tool name |
| `--length` | `8-20-4` | positive integer values |
| `--order` | `1-4-1` | used by BaMM |
| `--lpd` | `10-40-10` | used by SiteGA |
| `--nmotifs` | `5` | positive integer |
| `--jobs` | `1` | positive integer |
| `--seed` | unset | deterministic task seeds only when set |
| `--metric` | `pauROC` | `auPRC`, `auROC`, `pauPRC`, `pauROC` |
| `--comparison-criterion` | `score` | inclusive `>= 0.9` score or `<= 0.05` p-value by default |

P-value filtering requires the Mimosa null-distribution bundle and the MIMOSA
comparator. Nonzero external process status, missing output, invalid model
format, and incompatible model length are errors rather than empty success.

## Determinism

Odd/even indices have 1-based semantics. For task index `i`, a configured base
seed produces `base_seed + i + 1`; an unset seed is not replaced with a fake
deterministic value. Results are sorted by task index after parallel discovery.

## SiteGA Manifest

The required fields are:

```json
{
  "schema_version": 1,
  "status": "success",
  "message": "...",
  "models": [
    {
      "model_file": "model.mat",
      "model_type": "sitega",
      "length": 8,
      "name": "stable-name"
    }
  ]
}
```

MotifHORDE uses manifest order, validates the declared type and length against
the Mimosa model, and limits the returned vector to the requested count. It
does not infer semantics from filenames.
