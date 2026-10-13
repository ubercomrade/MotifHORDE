# Migration to MotifHORDE 0.3

MotifHORDE 0.3 is built on Mimosa 2.0 (`mimosa-tool>=2.0.0,<2.1`).

- Mimosa 2.0 has a single profile metric, weighted Jaccard. `--mimosa-metric`
  (`co`, `dice`, `cosine`) and the `mimosa_metric` argument of `DeNovoPipeline`
  and `ComparisonSettings` are removed; `compare_records` and
  `validate_null_distribution` no longer take `metric`.
- Motif comparison runs on a fixed, seeded set of 10,000 random sequences of
  length 100 (`comparison_sequences()`) instead of the foreground peaks. Scores
  no longer depend on the input data.
- Comparison scores are not comparable with 1.4 scores. Profiles keep only
  positions with `-log10(ERR) >= 4` (`min_logerr=4`, previously 0), and ERR is
  empirical on the compared sequences. Re-check `--comparison-threshold` (default still 0.9).
- Null bundles written by Mimosa 1.4 (format 9) are not compatible. Rebuild
  them with `mimosa.build_null` using `search_range` equal to
  `--mimosa-search-range`, `min_logerr=4`, and
  `sequences=comparison_sequences()`. The
  bundle contract is checked for metric `jaccard`, `alignment_version`,
  `search_range`, `min_logerr`, normalization, sequences, and background. The
  pipeline validates it once, in `DeNovoPipeline.__init__`.
- Model bundles (format 2) are unchanged and remain readable.

# Migration to MotifHORDE 0.2

MotifHORDE 0.2 is a breaking release built on Mimosa 1.4.

- `GenericModel`, the mutable model `config`, local scanning/model parsers,
  `TomtomComparator`, and pickle/joblib loaders are removed. Use Mimosa's
  public `MotifModel`, `scan`, `read_model`, and `write_model` APIs.
- Discovery and final models are Mimosa directory bundles, not `.pkl` files.
  They are stored under `output/<tool>/bootstrap/models` and
  `output/<tool>/motifs/models`.
- `--mimosa-null-distribution` accepts a Mimosa null-distribution directory
  bundle. Its metric, profile parameters, and normalization must match the
  comparison call; p-value BH correction is applied per query's target set.
- `compare_records` no longer revalidates the null bundle on every call. The
  pipeline validates the bundle once in `DeNovoPipeline.__init__` (without
  sequences) and once in `run()` (against the foreground sequences). Code
  that calls `compare_records` directly with a foreign bundle must run
  `validate_null_distribution` itself.
- MEME output is read directly by Mimosa. MotifHORDE no longer parses MEME
  matrices.
- SiteGA is an external process. It must write the versioned manifest passed by
  `--manifest`, including `schema_version`, `status`, `message`, and an ordered
  `models` list with `model_file`, `model_type`, `length`, and `name`.
- The CLI still takes four positional paths and supports `streme`, `meme`,
  `bamm`, `dimont`, `slim`, and `sitega`. Dimont/Slim JARs and all discovery
  programs remain external dependencies.

Regenerate old results from discovery inputs or interoperable MEME files.
