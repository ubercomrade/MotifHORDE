# CLI comparator and parallelism update plan

## Goal

Replace the old `continuous` comparator CLI with an explicit MIMOSA comparator
interface and make parallelism controlled by one shared CLI option.

Backward compatibility is not required. Do not keep deprecated aliases for:

- `-c continuous`;
- `--c-metric`;
- `--c-filter`;
- `--c-threshold`;
- `--c-search-range`;
- `--c-jobs`;
- `--tomtom-jobs`;
- `--meme-p`;
- `--jstacs-threads`.

The final interface should be smaller, explicit, and consistent across
discovery, comparison, and external tool execution.

## Current problems

### Comparator naming

The CLI exposes:

```text
-c, --comparator {tomtom,continuous}
```

The `continuous` name is stale. The implementation is already backed by
`mimosa.comparison`, and the user-facing comparator should be named `mimosa`.

### Comparison filtering

The current continuous comparator options expose:

```text
--c-filter {score,none}
--c-threshold FLOAT
```

This is too vague for the desired behavior. The comparison criterion should be
one of:

- `score`;
- `p-value`.

In the CLI, `p-value` means adjusted p-value. In MIMOSA result payloads the
public dataframe column is `adj.p-value`, not `p-value`.

### P-value configuration

MIMOSA p-values require a prepared null distribution. The CLI currently has no
way to provide that distribution to the pipeline comparator.

The existing `UniversalMotifComparator` constructor already accepts:

- `pvalue`;
- `null_distribution`;
- `null_search_dirs`;
- `effective_number_of_targets`.

However, the wrapper currently calls low-level `mimosa.comparison` functions.
Those low-level functions return score-only `ComparisonResult` records. The
p-value annotation logic is applied at the higher-level MIMOSA API boundary.

### Parallelism

Parallelism is currently split across multiple options:

- `--tomtom-jobs` for TomTom-like motif comparison;
- `--c-jobs` for continuous profile comparison;
- `--meme-p` for MEME;
- `--jstacs-threads` for Dimont and SlimDimont.

SiteGA has a `num_threads` binding parameter, but the CLI does not pass a thread
count into `SitegaDiscoveryTool`.

This creates inconsistent behavior and makes it unclear which option controls
which work.

## Target CLI

### Comparator selection

Replace comparator choices with:

```text
-c, --comparator {tomtom,mimosa}
```

Default remains:

```text
tomtom
```

### Shared parallelism

Add one shared option, preferably in `Other options`:

```text
--jobs INT
```

Suggested default:

```text
1
```

Semantics:

- `1` means sequential or single-threaded execution where supported;
- positive values mean that exact number of workers or threads where supported;
- `-1` means automatic all-core behavior.

Implementation detail:

- pass `args.jobs` directly to MIMOSA comparators, because MIMOSA normalizes
  `-1` to its automatic mode;
- for external tools that need a positive thread count, resolve `-1` to
  `os.cpu_count() or 1`;
- pass `None` only to APIs where omitting the thread option is intentional.

Add a small helper in `cli.py`:

```python
def _external_thread_count(jobs: int) -> int:
    if jobs == -1:
        return os.cpu_count() or 1
    return jobs
```

Validate `--jobs` once after parsing:

- valid values are `-1` and positive integers;
- reject `0` and values below `-1` with `parser.error(...)`.

### TomTom comparator options

Keep:

```text
--tomtom-metric {pcc,ed}
--pfm-mode
```

Remove:

```text
--tomtom-jobs
```

`TomtomComparator` should receive:

```python
n_jobs=args.jobs
```

### MIMOSA comparator options

Replace the old continuous group with a MIMOSA group.

Suggested options:

```text
--mimosa-metric {co,co_rowwise,dice,dice_rowwise,cosine}
--comparison-criterion {score,p-value}
--comparison-threshold FLOAT
--mimosa-search-range INT
--mimosa-null-distribution PATH
```

Optional MIMOSA null-distribution search support can be added if needed:

```text
--mimosa-null-search-dir PATH
--mimosa-effective-number-of-targets INT
```

Keep the first implementation minimal unless there is a current use case for
search directories or explicit E-value target count. The required parameter for
`p-value` mode is `--mimosa-null-distribution`.

Defaults:

- `--mimosa-metric co`;
- `--comparison-criterion score`;
- `--comparison-threshold` default should be derived from the criterion;
- `--mimosa-search-range 10`.

Threshold defaults:

- `score`: `0.9`;
- `p-value`: `0.05`.

Do not use one hard-coded threshold default in argparse, because the correct
default depends on the selected criterion. Set `default=None` in argparse and
resolve the final threshold after parsing.

Validation:

- if `--comparison-criterion p-value` is selected, require
  `--mimosa-null-distribution`;
- if `--comparison-criterion score` is selected, do not enable p-value
  annotation;
- reject missing `adj.p-value` result columns when p-value filtering is active.

## Code changes

### `src/motifhorde/cli.py`

Update parser examples:

- replace `-c continuous --c-metric co` with
  `-c mimosa --mimosa-metric co`;
- add one example for adjusted p-value filtering with
  `--comparison-criterion p-value` and `--mimosa-null-distribution`.

Update comparator choices:

```python
choices=["tomtom", "mimosa"]
```

Remove old comparator options:

- `--tomtom-jobs`;
- all `--c-*` options.

Add `--jobs` once.

Add a post-parse validation function, for example:

```python
def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    ...
```

Keep it focused on cross-option validation:

- `--jobs`;
- `--comparison-criterion p-value` requiring
  `--mimosa-null-distribution`;
- threshold default resolution.

Avoid scattering this validation through setup functions.

Update `setup_discovery_tool(args)`:

- MEME gets `threads=_external_thread_count(args.jobs)`;
- Dimont gets `threads=_external_thread_count(args.jobs)`;
- SlimDimont gets `threads=_external_thread_count(args.jobs)`;
- SiteGA gets `threads=_external_thread_count(args.jobs)`;
- STREME and BaMM remain unchanged unless a supported thread option is added
  deliberately.

Update `setup_comparator(args)`:

TomTom:

```python
return TomtomComparator(
    metric=args.tomtom_metric,
    n_jobs=args.jobs,
    seed=args.seed,
    pfm_mode=args.pfm_mode,
)
```

MIMOSA:

```python
return UniversalMotifComparator(
    name="mimosa_comparator",
    metric=args.mimosa_metric,
    n_jobs=args.jobs,
    seed=args.seed,
    comparison_criterion=args.comparison_criterion,
    comparison_threshold=args.comparison_threshold,
    search_range=args.mimosa_search_range,
    pvalue=args.comparison_criterion == "p-value",
    null_distribution=args.mimosa_null_distribution,
)
```

The exact parameter names on `UniversalMotifComparator` can be adjusted, but
avoid keeping `filter_type` as a public constructor concept if the pipeline now
uses a comparison criterion. `filter_type` is tied to the old `--c-filter`
model and is less explicit.

### `src/motifhorde/comparison.py`

Keep the wrappers thin. Do not duplicate MIMOSA null distribution loading or
p-value correction logic locally.

Change `UniversalMotifComparator.compare(...)` to use the public MIMOSA API
that applies p-value annotation:

```python
from mimosa.api import compare_one_to_many as mimosa_compare_one_to_many
```

Inside the wrapper, call the API with in-memory models:

```python
mimosa_compare_one_to_many(
    query=motif,
    targets=list(motifs_2),
    strategy="profile",
    sequences=sequences,
    comparator=self.config,
)
```

For `TomtomComparator`, evaluate whether it should also use the MIMOSA API with
`strategy="motif"`. This is preferable for consistency, but only change it if
tests confirm the output remains equivalent. If not, keep TomTom score-only for
this update and still pass `n_jobs=args.jobs`.

Normalize result frames through one helper:

```python
def _records_to_frame(records):
    return pd.DataFrame.from_records(records)
```

MIMOSA `ComparisonResult` is mapping-like and emits public columns such as:

- `score`;
- `p-value`;
- `adj.p-value`;
- `E-value`.

Do not rename MIMOSA output columns globally. Instead, map CLI criterion to the
column used by MotifHORDE filtering.

Suggested helper:

```python
def comparison_column_for_criterion(criterion: str) -> str:
    if criterion == "score":
        return "score"
    if criterion == "p-value":
        return "adj.p-value"
    raise ValueError(...)
```

### `src/motifhorde/pipeline.py`

The pipeline currently infers the comparison column from dataframe columns:

```python
def _comparison_column(frame: pd.DataFrame) -> str:
    if "p-value" in frame.columns:
        return "p-value"
    if "score" in frame.columns:
        return "score"
```

Replace implicit inference with explicit comparator configuration.

Add a small comparator-facing contract:

- comparator has `comparison_criterion`;
- comparator has `comparison_threshold`;
- comparator can expose `comparison_column`.

Keep this simple. Do not introduce abstract base classes or config objects for
this unless the code becomes materially clearer.

Filtering rules:

- `score`: similar if `score >= threshold`;
- `p-value`: similar if `adj.p-value <= threshold`.

Sorting rules:

- `score`: descending;
- `p-value`: ascending.

Update these functions to use the explicit column and direction:

- `_filter_similar_matches`;
- `_sort_comparisons`;
- `_is_similar_value`;
- `_comparison_sort_ascending`;
- `_comparison_column`.

Prefer passing the comparator or criterion into these helpers rather than using
module-level constants. The old constants:

```python
SIMILARITY_PVALUE_THRESHOLD = 0.001
SIMILARITY_SCORE_THRESHOLD = 0.9
```

should either be removed or replaced by named defaults in one place. Avoid
keeping stale constants that no longer control behavior.

Places that call `_filter_similar_matches` must be updated:

- `_compare_bootstrap_motifs`;
- `_select_best_full_motif`;
- `_select_nonredundant_motifs`;
- `_deduplicate_final_motifs`.

When p-value mode is active and `adj.p-value` is missing, raise a clear error:

```text
Comparison results do not contain adjusted p-values. Provide a compatible MIMOSA
null distribution with --mimosa-null-distribution.
```

### `src/motifhorde/discovery.py`

Update `SitegaDiscoveryTool`:

- add `threads: int | None = None` to `__init__`;
- store `self.threads`;
- pass `num_threads=self.threads or 0` to `sitega.train(...)`.

Use `0` only if the SiteGA binding contract still means "use OMP_NUM_THREADS".
If the CLI always passes a resolved positive value, passing `self.threads`
directly is simpler.

Do not add thread parameters to tools that do not support them.

Current supported mappings:

- MEME: `-p`;
- Dimont: `threads=...`;
- SlimDimont: `threads=...`;
- SiteGA: `num_threads=...`.

STREME and BaMM remain unchanged in this update.

## Documentation changes

Update `README.md`:

- replace `-c continuous` with `-c mimosa`;
- replace `--c-metric` with `--mimosa-metric`;
- document `--comparison-criterion score|p-value`;
- document that CLI `p-value` means MIMOSA `adj.p-value`;
- document that p-value mode requires `--mimosa-null-distribution`;
- document `--jobs` as the single parallelism option;
- remove references to `--tomtom-jobs`, `--c-jobs`, `--meme-p`, and
  `--jstacs-threads`.

Add a short example:

```bash
motifhorde peaks.fa bg.fa promoters.fa output/ \
  -c mimosa \
  --mimosa-metric co \
  --comparison-criterion p-value \
  --mimosa-null-distribution profile-null.joblib \
  --jobs -1
```

## Tests

### CLI parser tests

Update `tests/test_cli_tools.py` or add focused CLI tests:

- parser accepts `-c mimosa`;
- parser rejects `-c continuous`;
- help includes `--mimosa-metric`, `--comparison-criterion`,
  `--mimosa-null-distribution`, and `--jobs`;
- help does not include removed options;
- `--comparison-criterion p-value` without `--mimosa-null-distribution` fails;
- `--jobs 0` fails;
- `--jobs -2` fails.

### Comparator setup tests

Add tests for `setup_comparator(args)`:

- TomTom receives `n_jobs=args.jobs`;
- MIMOSA receives `n_jobs=args.jobs`;
- MIMOSA score mode sets `pvalue=False`;
- MIMOSA p-value mode sets `pvalue=True` and passes the null distribution path;
- MIMOSA p-value criterion maps to `adj.p-value`.

### Pipeline filtering tests

Add or update tests around the pure filtering helpers:

- score mode keeps rows with `score >= threshold`;
- score mode sorts descending;
- p-value mode keeps rows with `adj.p-value <= threshold`;
- p-value mode sorts ascending;
- p-value mode raises a clear error when `adj.p-value` is missing.

Keep these tests small and dataframe-based. Do not require real MIMOSA null
distribution files for pipeline filtering behavior.

### Discovery thread propagation tests

Update existing discovery tests:

- MEME command includes `-p <jobs>`;
- Dimont args include `threads=<jobs>`;
- SlimDimont args include `threads=<jobs>`;
- SiteGA `sitega.train(...)` receives `num_threads=<jobs>`.

Add CLI setup tests:

- `--jobs 3 -t meme` creates `MemeDiscoveryTool` with `threads == 3`;
- `--jobs 3 -t dimont` creates `DimontDiscoveryTool` with `threads == 3`;
- `--jobs 3 -t slim` creates `SlimDiscoveryTool` with `threads == 3`;
- `--jobs 3 -t sitega` creates `SitegaDiscoveryTool` with `threads == 3`.

For `--jobs -1`, monkeypatch `os.cpu_count()` and assert external tool thread
counts receive the resolved positive value.

### MIMOSA p-value integration test

Add one focused test around `UniversalMotifComparator` with monkeypatched MIMOSA
API:

- monkeypatch `motifhorde.comparison.mimosa_compare_one_to_many`;
- return one `ComparisonResult` or dict containing `adj.p-value`;
- assert `UniversalMotifComparator.compare(...)` returns a frame with
  `adj.p-value`.

This verifies MotifHORDE uses the API boundary that can annotate p-values
without requiring an expensive real null distribution fixture.

## Implementation order

1. Update CLI parser and argument validation.
2. Add shared `--jobs` plumbing into comparator setup.
3. Add shared `--jobs` plumbing into discovery setup.
4. Extend `SitegaDiscoveryTool` to accept and pass thread count.
5. Update `UniversalMotifComparator` to use MIMOSA API for one-to-many
   comparison.
6. Replace implicit pipeline comparison-column detection with explicit criterion
   and threshold handling.
7. Update tests for parser, setup, filtering, and thread propagation.
8. Update README examples and option documentation.
9. Run formatting and tests.

Suggested verification commands:

```bash
uv run ruff check .
uv run pytest
```

If external smoke tests are needed, run them separately because they depend on
installed external tools and are slower.

## Expected behavior changes

Breaking CLI changes:

- `-c continuous` is invalid;
- all `--c-*` options are invalid;
- `--tomtom-jobs` is invalid;
- `--meme-p` is invalid;
- `--jstacs-threads` is invalid.

New behavior:

- `-c mimosa` selects profile comparison through MIMOSA;
- `--jobs` controls all supported comparison and discovery parallelism;
- score filtering uses `score >= threshold`;
- p-value filtering uses `adj.p-value <= threshold`;
- p-value mode requires a prepared MIMOSA null distribution.

## Compatibility risks

The main risk is that existing scripts using old CLI options will fail. This is
acceptable because backward compatibility is explicitly not required.

The second risk is p-value support through MIMOSA API. This should be tested
with a monkeypatched unit test first, then with one real compatible null
distribution file if a fixture is available.

The third risk is `--jobs -1` behavior for external tools. Resolve it once in
CLI setup and pass only positive thread counts to tools that require positive
values.
