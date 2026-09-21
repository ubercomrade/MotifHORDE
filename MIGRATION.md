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
