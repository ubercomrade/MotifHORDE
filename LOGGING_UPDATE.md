# Logging Update Plan

## Goal

Make CLI stdout more informative and easier to inspect during real pipeline runs.
The output should clearly show:

- what stage is running;
- which discovery tool, comparator, metric, and parameters are active;
- how many inputs, tasks, motifs, comparisons, matches, and filtered records were
  produced;
- where output files are written;
- how long important stages took.

The implementation should stay simple. Use the Python standard library
`logging` module, avoid a custom logging framework, and keep pure computation
helpers free of side effects.

## Target stdout format

Use one event per line:

```text
2026-06-26 18:42:11 INFO  [pipeline] start | tool=sitega | comparator=tomtom | metric=pauROC | jobs=1 | seed=42
2026-06-26 18:42:11 INFO  [input] loaded | foreground=400 | background=400 | promoters=400
2026-06-26 18:42:11 INFO  [bootstrap] start | param_sets=1 | tasks=2 | motifs_per_task=1 | fpr=0.001
2026-06-26 18:42:18 INFO  [bootstrap.discovery] done | task=1/2 | split=odd | params=length:6,lpd:10 | motifs_found=1 | elapsed=7.1s
2026-06-26 18:42:26 INFO  [bootstrap.compare] done | params=length:6,lpd:10 | odd=1 | even=1 | odd_selected=1 | even_selected=1 | raw_pairs=1 | passed=1 | deduped=1
2026-06-26 18:43:01 INFO  [final.discovery] done | params=length:6,lpd:10 | requested_motifs=2 | motifs_found=2 | elapsed=35.2s
2026-06-26 18:43:02 INFO  [final.match] done | params=length:6,lpd:10 | candidate_pairs=1 | matched=1 | unmatched=0
2026-06-26 18:43:02 INFO  [final.dedup] done | candidates=1 | kept=1 | removed=0
2026-06-26 18:43:02 INFO  [output] saved | bootstrap_models=2 | final_models=1 | dir=/path/out/sitega
2026-06-26 18:43:02 INFO  [pipeline] done | elapsed=51.3s
```

Format rules:

- Timestamp: `YYYY-MM-DD HH:MM:SS`.
- Level: left-aligned, for example `INFO ` or `WARN `.
- Stage tag: short bracketed name such as `[bootstrap.compare]`.
- Message: `event | key=value | key=value`.
- Parameters: stable compact format, for example `params=length:6,lpd:10`.
- Durations: seconds with one decimal place, for example `elapsed=7.1s`.
- Avoid nested, indented multiline status blocks in normal verbose output.

## Design Decisions

1. Use stdlib `logging`.

   This avoids custom stdout plumbing and lets existing `logger.warning(...)`
   calls in `discovery.py` participate in the same output format.

2. Configure logging once in `cli.py`.

   CLI owns presentation. Core code should only emit structured events through
   module loggers.

3. Keep verbose behavior explicit.

   `--verbose` should set the log level to `INFO`. Without `--verbose`, normal
   runs should stay quiet except warnings and errors.

4. Keep computation helpers pure.

   Do not add logging inside helpers such as `_filter_similar_matches`,
   `_deduplicate_matches`, `_select_nonredundant_motifs`, or
   `_deduplicate_final_motifs` unless the function is changed to return
   explicit counters. Prefer logging at the orchestration call sites where
   stage context is available.

5. Do not log from bootstrap worker processes.

   Bootstrap discovery can use `ProcessPoolExecutor`. Logging from worker
   processes may interleave output and make logs harder to read. Workers should
   return data; the parent process should log completed task summaries.

## Implementation Steps

### 1. Add logging setup in `cli.py`

Add a small setup function near CLI utilities:

```python
def configure_logging(verbose: bool) -> None:
    level = logging.INFO if verbose else logging.WARNING
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)-5s [%(name)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stdout,
        force=True,
    )
```

Call it immediately after argument validation in `main_cli()`:

```python
args = parser.parse_args()
validate_args(parser, args)
configure_logging(args.verbose)
```

Keep user-facing fatal CLI validation errors as direct `print()` plus `sys.exit`
for now, because they happen before the pipeline starts and are simple.

### 2. Replace verbose CLI header with log events

Replace the current block in `main_cli()` that prints:

- `MotifHORDE De Novo Pipeline`;
- discovery tool;
- comparator;
- metric;
- jobs;
- bootstrap workers.

Use one or two `logging.getLogger("pipeline").info(...)` calls instead:

```text
[pipeline] start | tool=... | comparator=... | metric=... | fpr=... | jobs=...
[pipeline] workers | bootstrap_discovery=... | discovery_threads=1 | comparator_jobs=...
```

This keeps the startup information but removes separator lines and manual
formatting.

### 3. Add small formatting helpers

Add helpers in `pipeline.py` or a tiny new module only if reuse justifies it.
The simplest first pass is to keep them in `pipeline.py`:

```python
def _format_log_params(params: dict[str, Any]) -> str:
    return ",".join(f"{key}:{params[key]}" for key in sorted(params))


def _format_elapsed(seconds: float) -> str:
    return f"{seconds:.1f}s"
```

Do not introduce classes for formatting. These are small pure functions with
clear reuse.

### 4. Add pipeline-level logger

In `pipeline.py`:

```python
import logging
import time

logger = logging.getLogger("pipeline")
```

Use child loggers only when they improve tags:

```python
input_logger = logging.getLogger("input")
bootstrap_logger = logging.getLogger("bootstrap")
compare_logger = logging.getLogger("bootstrap.compare")
final_logger = logging.getLogger("final.discovery")
```

Alternative: keep one module logger and include the stage in the message. The
preferred project format uses logger names as tags, so separate named loggers
are clearer.

### 5. Log input loading

In `DeNovoPipeline.run()`, after `_read_sequences(...)`:

```text
[input] loaded | foreground=... | background=... | promoters=...
```

Use `len(batch["lengths"])` for counts.

Also log prepared output directories after `_prepare_output_dirs(...)`:

```text
[output] prepared | bootstrap_dir=... | motifs_dir=...
```

### 6. Pass logging context to `Bootstrapper`

Add `verbose` only if needed for compatibility, but prefer relying on global
logging level instead of storing `verbose` in every class.

Change `Bootstrapper.run()` to emit `INFO` logs. When logging is configured at
`WARNING`, those messages are automatically suppressed.

No `verbose` attribute is required if all output goes through logging.

### 7. Instrument bootstrap task construction and discovery

In `Bootstrapper.run()`:

- record `start_time = time.perf_counter()`;
- call `_build_bootstrap_tasks(...)`;
- log number of parameter sets, tasks, workers, requested motifs, FPR threshold;
- run discovery tasks;
- log one summary line per completed result in parent process;
- log total discovery elapsed.

Recommended messages:

```text
[bootstrap] start | param_sets=3 | tasks=6 | jobs=2 | motifs_per_task=5 | fpr=0.001
[bootstrap.discovery] done | task=1/6 | split=odd | params=length:8 | motifs_found=5
[bootstrap.discovery] done | task=2/6 | split=even | params=length:8 | motifs_found=5
[bootstrap.discovery] complete | tasks=6 | motifs_found=28 | elapsed=83.4s
```

For elapsed per task in parallel mode, either:

- skip per-task elapsed in the first implementation; or
- add `elapsed` to `BootstrapDiscoveryResult` from inside `_run_discovery_task`.

Adding `elapsed` to the returned result is acceptable because it is explicit
data, not worker stdout.

### 8. Instrument bootstrap evaluation

In `_evaluate_bootstrap_results(...)`, count evaluated motifs. Since this
function currently returns only statistics and motifs, keep it side-effect-free
or log only in `Bootstrapper.run()` after receiving the return value.

Preferred simple approach:

```python
statistics, bootstrap_motifs = _evaluate_bootstrap_results(...)
logging.getLogger("bootstrap.evaluate").info(
    "done | motifs=%d | statistics=%d",
    len(bootstrap_motifs),
    len(statistics),
)
return statistics, bootstrap_motifs
```

This keeps evaluation logging in the orchestration method.

### 9. Replace bootstrap save `print()`

In `_save_bootstrap(...)`, replace:

```python
print(f"Saving {len(motifs)} bootstrap motifs to {models_dir}...")
```

with:

```text
[output] saved | type=bootstrap_models | count=... | dir=...
[output] saved | type=bootstrap_statistics | path=...
```

This method is I/O shell code, so logging there is appropriate.

### 10. Instrument bootstrap motif comparison

In `_compare_bootstrap_motifs(...)`, split the current pipeline into named
intermediate frames:

```python
raw_frame = self.comparator.compare(...)
similar_frame = _filter_similar_matches(raw_frame, self.comparator)
sorted_frame = _sort_comparisons(similar_frame, self.comparator)
deduped_frame = _deduplicate_matches(sorted_frame)
```

Log per parameter set:

```text
[bootstrap.compare] start | params=length:8 | odd=5 | even=5
[bootstrap.compare] selected | params=length:8 | odd_selected=3 | even_selected=4
[bootstrap.compare] done | params=length:8 | raw_pairs=12 | passed=5 | deduped=3
```

If a parameter set is skipped, log the reason:

```text
[bootstrap.compare] skipped | params=length:8 | reason=missing_odd_or_even | odd=0 | even=5
[bootstrap.compare] skipped | params=length:8 | reason=no_nonredundant_motifs | odd_selected=0 | even_selected=2
[bootstrap.compare] skipped | params=length:8 | reason=no_pairs_passing_threshold | raw_pairs=12 | passed=0
```

At the end:

```text
[bootstrap.compare] complete | param_sets=3 | matched_param_sets=2 | records=7
```

### 11. Instrument final discovery

In `_select_final_motifs(...)`, around `self.discovery_tool.discover(...)`:

```text
[final.discovery] start | params=length:8 | requested_motifs=10
[final.discovery] done | params=length:8 | motifs_found=8 | elapsed=35.2s
```

Also log how many bootstrap comparison records will be used:

```text
[final.match] start | params=length:8 | candidate_pairs=3
```

Track counters in the existing loop:

- `matched`;
- `unmatched`;
- `already_assigned` if a candidate is skipped because no full motifs remain.

Log:

```text
[final.match] selected | params=length:8 | pair=odd_name/even_name | motif=Full-1
[final.match] unmatched | params=length:8 | pair=odd_name/even_name
[final.match] done | params=length:8 | candidate_pairs=3 | matched=2 | unmatched=1
```

The per-pair `selected` line can be useful but may be noisy. If logs become too
long, keep only the final `done` line and the existing final motif summary.

### 12. Instrument best full motif selection only through returned data

Do not put logging inside `_select_best_full_motif(...)` initially. It is mostly
selection logic and is easier to test without side effects.

If detailed diagnostics are needed later, change it to return a small
`TypedDict` with:

- selected motif;
- compared candidates;
- candidates passing odd reference;
- candidates passing even reference;
- final candidates.

Do not add hidden mutation or global diagnostic state.

### 13. Instrument final deduplication

Before and after `_deduplicate_final_motifs(...)` in `run()`:

```text
[final.dedup] start | candidates=...
[final.dedup] done | candidates=... | kept=... | removed=...
```

Do not log inside `_deduplicate_final_motifs(...)` in the first pass.

### 14. Replace final save `print()`

In `_save_results(...)`, replace:

- `Saving ... individual models`;
- final motif summary prints.

Use:

```text
[output] saved | type=final_models | count=... | dir=...
[output] saved | type=final_meme | path=...
[output] saved | type=final_statistics | path=...
[motif] rank=1 | name=... | params=length:8 | auPRC=... | auROC=... | pauPRC=... | pauROC=...
```

Use logger name `motif` for ranked motif summaries.

### 15. Normalize warnings

`discovery.py` already uses `logger.warning(...)` for invalid or missing
outputs. After `configure_logging(...)`, these warnings will share the same
timestamp/level/tag format.

Replace the direct warning-like print in `PerformanceEvaluator.evaluate(...)`:

```python
print(f"Incorrect background_type: {self.background_type}, set as `peaks`")
```

with:

```python
logger.warning(
    "invalid_background_type | value=%s | fallback=peaks",
    self.background_type,
)
```

However, because CLI currently restricts background type choices, this is mostly
a defensive path.

### 16. Keep external command stdout/stderr behavior unchanged

`run_checked(...)` captures external stdout/stderr and includes it in raised
errors. Do not stream external tool output into normal verbose logs in the first
implementation.

Reason:

- external tools can produce large, inconsistent output;
- current fullrun logs already capture pipeline stdout/stderr;
- adding external streaming would make normal INFO logs harder to read.

If needed later, add a separate `--debug-external` option or write external
stdout/stderr to files under the task output directory.

## Testing Plan

### Unit tests

Add tests for formatting helpers:

- `_format_log_params({"length": 6, "lpd": 10}) == "length:6,lpd:10"`;
- stable sorted key order;
- elapsed formatting.

### CLI logging tests

Add or update tests that run the CLI with `--verbose` on a small fake or existing
small dataset and assert stdout contains:

- `INFO`;
- `[pipeline] start`;
- `[input] loaded`;
- `[bootstrap] start`;
- `[bootstrap.compare]`;
- `[final.dedup]`;
- `[pipeline] done`.

Also assert no old separator-heavy header is required by tests except the
external fullrun test. Update that test to look for `[pipeline] start` instead
of `MotifHORDE De Novo Pipeline`.

### Non-verbose tests

Run a small CLI path without `--verbose` and assert normal INFO events are not
printed. Warnings may still appear.

### Bootstrap tests

Existing bootstrap tests should continue to pass. If `BootstrapDiscoveryResult`
gets a new `elapsed` key, update expected structures only if tests compare full
dicts. Current tests mostly compare returned motif names and stats, so impact
should be low.

### Full test commands

Fast validation:

```bash
uv run pytest tests/test_bootstrap_parallel.py tests/test_pipeline_selection.py tests/test_cli_tools.py -q
```

Broader validation:

```bash
uv run pytest -q
```

External real-data validation remains opt-in and should only be run when needed:

```bash
HORDEMOTIFS_RUN_FULLRUN=1 \
HORDEMOTIFS_FULLRUN_LOG_DIR=/tmp/hordemotifs-real-fullrun-logs \
HORDEMOTIFS_STREME_COMMAND=/home/anton/miniconda3/envs/motifhorde/bin/streme \
HORDEMOTIFS_JAVA_COMMAND=/home/anton/miniconda3/bin/java \
/home/anton/miniconda3/envs/motifhorde/bin/python -m pytest \
  tests/test_external_fullrun.py::test_real_data_full_pipeline_smoke_with_verbose_log \
  -q
```

## Rollout Order

1. Add `configure_logging(...)` in `cli.py`.
2. Add formatting helpers and module loggers.
3. Convert CLI verbose header/footer to log events.
4. Convert unconditional `print()` calls in `pipeline.py` and `evaluation.py`.
5. Add bootstrap discovery/evaluation counters.
6. Add bootstrap comparison counters.
7. Add final discovery/matching/dedup/output counters.
8. Update tests that assert old stdout text.
9. Run fast tests.
10. Optionally run external fullrun smoke test and inspect log readability.

## Compatibility Notes

- `--verbose` remains the user-facing switch for detailed output.
- Without `--verbose`, INFO logs are suppressed.
- Warnings from `discovery.py` become more visible and consistently formatted.
- Fullrun logs will change because the old `MotifHORDE De Novo Pipeline` banner
  should be replaced by `[pipeline] start`.
- Existing programmatic API users may see warnings if they configure root logging
  themselves, but pipeline INFO output remains controlled by logging level.

## Risks

- Switching from `print()` to `logging` may require updating tests that inspect
  stdout.
- `logging.basicConfig(force=True)` can override application logging when
  `main_cli()` is called from another Python process. This is acceptable for a
  CLI entry point, but should not be used from library code.
- Adding too many per-pair logs could make real-data fullrun logs noisy. Prefer
  summary counters first, then add detailed per-pair logs only if debugging
  requires them.

## Definition of Done

- `--verbose` output has timestamps, levels, stage tags, and key-value fields.
- Normal non-verbose runs do not print INFO progress events.
- Bootstrap logs show task counts, parameter sets, motifs found, and evaluation
  totals.
- Comparison logs show raw pairs, threshold-passing pairs, and deduplicated
  pairs.
- Final logs show discovery counts, matching counts, deduplication counts, and
  saved output paths.
- Direct `print()` calls in pipeline/evaluation progress paths are removed.
- Tests cover formatting and at least one verbose CLI stdout path.
- Fast test suite passes.
