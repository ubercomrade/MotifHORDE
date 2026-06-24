# Plan: bootstrap parallelism through the existing `--jobs`

## Goal

Speed up SiteGA-heavy runs without adding another jobs option and without
reintroducing unsafe OpenMP mutation inside the SiteGA C++ core.

Target example:

```bash
motifhorde peaks.fa background.fa promoters.fa output \
  --tool sitega \
  -n 10 \
  --lpd 10-40-10 \
  -l 10-16-2 \
  -m auROC \
  -c mimosa \
  --mimosa-metric dice \
  --jobs 6 \
  -v
```

For this grid the bootstrap stage currently runs at least:

```text
4 motif lengths * 4 lpd values * 2 odd/even splits = 32 discovery runs
```

Those discovery runs are independent. They are a safer and simpler
parallelization boundary than the shared-population mutation loop inside
SiteGA.

## User-facing decision

Do not add `--bootstrap-jobs`.

Use the existing `--jobs` as the single parallelism budget:

- during bootstrap, `--jobs N` means up to `N` independent discovery runs in
  separate Python processes;
- discovery tools must not use their own multithreading while bootstrap process
  parallelism is active;
- after bootstrap, motif comparison uses the same `--jobs N`;
- MIMOSA should use its existing internal `n_jobs` support;
- Tomtom should not get a new parallelism mechanism.

This keeps the CLI small and avoids asking users to reason about multiple job
counts.

## Current state

- `Bootstrapper.run()` in `src/motifhorde/evaluation.py` iterates over the
  parameter grid and odd/even splits sequentially.
- `setup_discovery_tool()` currently passes `_external_thread_count(args.jobs)`
  into tools that support thread options.
- `setup_comparator()` currently passes `args.jobs` into MIMOSA/Tomtom
  comparator configuration.
- SiteGA training is now serial inside the C++ core because the previous OpenMP
  mutation loop read and mutated shared `pop[]` state and could corrupt memory.
- SiteGA still accepts `num_threads`, but the current safe C++ path does not use
  OpenMP for training work.
- SiteGA uses `time(NULL)` as the default seed when no seed is passed. Parallel
  process launches can start in the same second, so explicit per-task seeds are
  required for reproducible and diverse SiteGA bootstrap runs.

## Target behavior

- `--jobs 1` preserves current sequential bootstrap behavior.
- `--jobs N`, where `N > 1`, runs independent bootstrap discovery tasks in a
  `ProcessPoolExecutor(max_workers=N)`.
- `--jobs -1` resolves to `os.cpu_count() or 1` and uses that value for
  bootstrap process workers and comparator jobs.
- Discovery tools run with one internal thread during bootstrap. This prevents
  `N` bootstrap workers from each starting `N` tool threads.
- MIMOSA comparison runs after bootstrap and keeps using its existing internal
  `n_jobs=args.jobs` path.
- Tomtom comparison is left unchanged. No new Tomtom-level process pool or
  extra parallelism should be added.
- Output ordering remains deterministic by sorting bootstrap results by task
  index.
- Default command behavior remains understandable: one option, one compute
  budget.

## Non-goals for the first patch

- Do not restore OpenMP inside SiteGA mutation.
- Do not rewrite the SiteGA genetic algorithm.
- Do not add `--bootstrap-jobs`, `--discovery-jobs`, or any other jobs option.
- Do not add a separate MIMOSA parallelism layer.
- Do not add a separate Tomtom parallelism layer.
- Do not parallelize final motif selection in the same patch.
- Do not introduce a generic executor framework or service layer.

## CLI changes

Keep the existing `--jobs` option and validation:

```text
--jobs INT
```

Semantics after this change:

- `1`: sequential bootstrap and comparator/internal jobs of 1;
- positive integer: bootstrap process worker count and comparator job count;
- `-1`: resolve to all available CPUs for bootstrap worker count. Comparator
  setup should keep its current handling of `args.jobs` unless a specific
  comparator already resolves `-1` internally.

Do not add another public argument.

Update help text to make the dual-stage behavior explicit:

```text
Shared worker count. During bootstrap, independent discovery runs use this many
processes and discovery tools run single-threaded inside each process. During
comparison, supported comparators use this many internal jobs. Use -1 for all
available cores.
```

Keep validation simple:

```python
if args.jobs == 0 or args.jobs < -1:
    parser.error("--jobs must be -1 or a positive integer")
```

Keep a single resolver:

```python
def _resolve_jobs(jobs: int) -> int:
    if jobs == -1:
        return os.cpu_count() or 1
    return jobs
```

Use the resolved value for bootstrap worker processes. Keep comparator setup on
the existing code path so Tomtom behavior does not change.

Verbose output should include:

```text
Jobs: 6
Bootstrap discovery workers: 6
Discovery tool internal threads during bootstrap: 1
Comparator jobs: 6
```

The exact wording can be shorter, but it should make the policy visible.

## Parallelism policy

Use `--jobs` as a global compute budget, not as a multiplier.

### Bootstrap discovery

- Use process-level parallelism.
- Worker count is the resolved `--jobs` value.
- Inside each bootstrap worker, discovery tools must use a single internal
  thread or no explicit thread option.
- This applies to tools with internal parallelism support: MEME, DiMotif,
  SlimDimont, SiteGA, and any future discovery tool that exposes thread count.

### Motif comparison

- Comparison happens after bootstrap.
- It can use the same resolved `--jobs` value because bootstrap workers are no
  longer running.
- MIMOSA should use its existing internal `n_jobs` implementation.
- Tomtom should remain as currently implemented. Do not add a new process pool
  around Tomtom.

### Final discovery

- Keep final discovery sequential in the first patch.
- Use discovery tools with one internal thread for consistency with the new
  discovery policy unless there is a strong reason to preserve internal threads
  for final discovery.
- If final discovery becomes a bottleneck later, handle it as a separate phase.

## Discovery tool construction

The current `setup_discovery_tool(args)` creates tools with
`threads=_external_thread_count(args.jobs)` for tools that support threads.
That must change to avoid oversubscription.

Target behavior:

- `setup_discovery_tool(args)` should configure discovery tools with internal
  thread count `1` or `None`, depending on each tool's existing API.
- The resolved `jobs` value should be passed to the pipeline/bootstrap layer,
  not into discovery tool internals.
- `setup_comparator(args)` should keep its existing behavior. In particular,
  MIMOSA should use its internal `n_jobs` path, and Tomtom should not receive a
  new wrapper-level parallelism mechanism.

Recommended implementation:

```python
def setup_discovery_tool(args) -> Any:
    discovery_threads = 1
    ...
    return MemeDiscoveryTool(..., threads=discovery_threads)
    ...
    return SitegaDiscoveryTool(..., threads=discovery_threads)
```

For tools where `None` means "tool default may use multiple cores", prefer `1`
if the tool supports an explicit single-thread option. Use `None` only for tools
that do not expose thread control.

Document this as an intentional policy, not an accidental limitation.

## Pipeline wiring

Add `jobs` to `DeNovoPipeline`:

```python
class DeNovoPipeline:
    def __init__(
        self,
        discovery_tool: MotifDiscoveryTool,
        evaluator: PerformanceEvaluator,
        comparator: UniversalMotifComparator,
        fpr_threshold: float = 0.001,
        number_of_motifs: int = 5,
        jobs: int = 1,
        seed: int | None = None,
    ) -> None:
        ...
```

Pass `jobs` and `seed` into `Bootstrapper`:

```python
bootstrapper = Bootstrapper(
    self.discovery_tool,
    self.evaluator,
    output_dir,
    jobs=self.jobs,
    seed=self.seed,
)
```

Use the existing CLI `--seed` as the base seed for bootstrap task seed
derivation. This keeps one reproducibility control for the whole pipeline.

Do not pass comparator objects into bootstrap workers. Comparison happens after
bootstrap in the parent process. The `jobs` value here is the resolved positive
worker count for bootstrap, not a new CLI option.

## Bootstrap task model

Use `TypedDict` instead of new data classes. This keeps the implementation
simple and aligned with the repository style.

Suggested shapes in `src/motifhorde/evaluation.py`:

```python
class BootstrapTask(TypedDict):
    index: int
    params: dict[str, Any]
    params_suffix: str
    step_name: str
    fg_path: str
    bg_path: str
    output_dir: str
    number_of_motifs: int
    seed: int | None


class BootstrapDiscoveryResult(TypedDict):
    index: int
    params: dict[str, Any]
    params_suffix: str
    step_name: str
    motifs: list[GenericModel]
```

Keep task data explicit. Do not hide paths, params, jobs, or seed in globals.

## Bootstrap implementation

Refactor `Bootstrapper.run()` into small functions with clear responsibility:

```text
run()
  prepare task root directory
  build bootstrap discovery tasks and write per-task FASTA files
  run discovery tasks sequentially or in a process pool
  evaluate returned motifs in the parent process
  return statistics and bootstrap motifs
```

Suggested private functions:

```python
def _bootstrap_indices(n_peaks: int, step_name: str) -> tuple[list[int], list[int]]:
    ...


def _build_bootstrap_tasks(...) -> tuple[list[BootstrapTask], dict[int, SequenceBatch]]:
    ...


def _run_discovery_task(
    discovery_tool: MotifDiscoveryTool,
    task: BootstrapTask,
) -> BootstrapDiscoveryResult:
    ...


def _run_bootstrap_discovery_tasks(
    discovery_tool: MotifDiscoveryTool,
    tasks: list[BootstrapTask],
    jobs: int,
) -> list[BootstrapDiscoveryResult]:
    ...


def _evaluate_bootstrap_results(...) -> tuple[dict[str, Any], list[GenericModel]]:
    ...
```

`_run_discovery_task` must be a top-level function, not a nested function, so it
can be used by `ProcessPoolExecutor` with spawn-style process creation.

The worker should do only discovery:

```python
def _run_discovery_task(discovery_tool, task):
    try:
        kwargs = dict(task["params"])
        if task["seed"] is not None and discovery_tool.name == "sitega":
            kwargs["seed"] = task["seed"]

        motifs = discovery_tool.discover(
            task["fg_path"],
            task["bg_path"],
            task["output_dir"],
            number_of_motifs=task["number_of_motifs"],
            **kwargs,
        )
    except Exception as exc:
        raise RuntimeError(
            f"Bootstrap discovery failed for params={task['params']} "
            f"split={task['step_name']}"
        ) from exc

    return {
        "index": task["index"],
        "params": task["params"],
        "params_suffix": task["params_suffix"],
        "step_name": task["step_name"],
        "motifs": motifs,
    }
```

Evaluate motifs in the parent process after discovery completes. This avoids
pickling and repeatedly sending the full test/background sequence batches to
workers. It also keeps side effects simpler:

```python
for result in sorted(results, key=lambda item: item["index"]):
    test_peaks = test_batches[result["index"]]
    for motif in result["motifs"]:
        stats = self.evaluator.evaluate(
            motif,
            test_peaks,
            background,
            err_threshold,
        )
        motif.name = (
            f"{motif.name}_{result['params_suffix']}_{result['step_name']}"
        )
        statistics[motif.name] = stats
        bootstrap_motifs.append(motif)
```

This keeps output names and result ordering deterministic.

## Process pool details

Use `ProcessPoolExecutor` only when `jobs > 1`:

```python
if jobs == 1:
    return [_run_discovery_task(discovery_tool, task) for task in tasks]

mp_context = multiprocessing.get_context("spawn")
with ProcessPoolExecutor(max_workers=jobs, mp_context=mp_context) as executor:
    results = list(
        executor.map(
            _run_discovery_task_with_tool,
            repeat(discovery_tool),
            tasks,
        )
    )
```

Implementation can use either:

- a small top-level wrapper accepting `(discovery_tool, task)`;
- `executor.submit(...)` with explicit arguments.

Prefer the simpler version that is easiest to test.

Use `spawn` for the parallel path. It is safer for native extensions and global
C++ state than inheriting already-imported state with `fork`. The startup
overhead is small relative to SiteGA training time.

## Temporary directory handling

Create one parent temporary directory for all bootstrap tasks:

```python
with tempfile.TemporaryDirectory(
    dir=os.path.join(self.output_dir, self.discovery_tool.name),
    prefix="bootstrap_parallel_",
) as task_root:
    ...
```

Inside it create stable per-task directories:

```text
task_0000_length-10_lpd-10_odd/
task_0001_length-10_lpd-10_even/
...
```

The parent process should write:

- `train.fasta`;
- `background.fasta`.

Each worker uses its own task directory as discovery output. This avoids file
name collisions from external tools and from SiteGA log/model files.

The temporary root must stay alive until:

- all workers finish;
- all returned motifs are evaluated;
- no worker can still read input FASTA files.

## Seed handling

Add explicit SiteGA seed support:

1. Add `seed: int | None = None` to `SitegaDiscoveryTool.__init__`.
2. Store it as `self.seed`.
3. In `discover()`, pass an explicit seed to `sitega.train`:

```python
seed = int(kwargs.get("seed", self.seed or 0))
...
sitega.train(..., seed=seed, ...)
```

Derive per-task seeds in `Bootstrapper` when a base seed is available:

```python
def _task_seed(base_seed: int | None, task_index: int) -> int | None:
    if base_seed is None:
        return None
    return base_seed + task_index + 1
```

Keep this deliberately simple and deterministic. Avoid Python's built-in
`hash()`, because it is salted per process.

If no base seed is provided, leave `seed=None` and keep current non-deterministic
SiteGA behavior.

Do not inject `seed` into every discovery tool in the first patch. Some tools
already manage seeds through constructor arguments, and some ignore seeds. For
the first patch, pass per-task seed only to SiteGA, where the C++ binding already
supports it and where parallel same-second `time(NULL)` collisions matter.

## Tool-specific impact

### SiteGA

- Main beneficiary.
- Each training run remains serial and memory-isolated.
- Multiple trainings can run concurrently through separate Python processes.
- `--jobs N` controls how many SiteGA trainings run at the same time during
  bootstrap.
- SiteGA receives `num_threads=1` in discovery tool setup.
- Recommended usage:

```bash
--jobs 6 --seed 42
```

### MEME and STREME

- Parallel bootstrap can run multiple external tool processes at once.
- Internal thread options should be disabled or set to `1` for discovery.
- This may improve wall time if CPU and RAM are available.
- Watch for external tool memory use and temporary output volume.

### BaMM

- Can benefit if individual BaMM runs are independent and system resources are
  sufficient.
- Avoid internal multithreading in bootstrap subprocesses where possible.

### DiMotif/SlimDimont

- Parallel bootstrap may launch multiple JVMs.
- Each JVM should be configured for one internal thread where the tool supports
  it.
- This can be memory-heavy because each JVM may reserve up to `java_xmx`.
- Users should choose `--jobs` conservatively for Java tools.

### MIMOSA comparison

- MIMOSA comparison happens after bootstrap, so it can use the same resolved
  `--jobs` budget.
- Do not add an outer MIMOSA process pool.
- Keep using the existing `UniversalMotifComparator(..., n_jobs=..., ...)` setup
  path and rely on MIMOSA's internal parallel implementation.

### Tomtom comparison

- Do not change Tomtom parallelism in this patch.
- Keep existing `TomtomComparator(..., n_jobs=args.jobs, ...)` behavior if that
  is already how Tomtom is configured.
- Do not add a new process pool or task scheduler around Tomtom.

## Final discovery parallelism

Do not parallelize final discovery in the first patch.

Reason:

- final selection depends on bootstrap comparison records;
- each final discovery is followed by matching against odd/even references;
- the logic is more coupled than bootstrap discovery;
- bootstrap is the larger and simpler win for SiteGA grids.

Possible phase 2:

- collect final parameter groups;
- run full discovery for each group in a process pool using the same `--jobs`
  budget;
- keep discovery tool internal threads at `1`;
- return full motif lists;
- perform matching and final motif assignment in the parent process.

Keep this as a separate patch after bootstrap parallelism is stable.

## Optional SiteGA speed knobs

These should be separate follow-up changes, not part of bootstrap parallelism:

### `--sitega-pop-size`

Expose the existing C++ `pop_size` parameter:

```text
--sitega-pop-size INT
```

Expected effect:

- lower values should reduce runtime;
- lower values may reduce motif quality or search stability.

Implementation:

- add `pop_size: int | None = None` to `SitegaDiscoveryTool`;
- pass `pop_size=self.pop_size or 0` to `sitega.train`;
- validate `1 <= pop_size <= 500` if provided.

### `--sitega-max-peak-len`

Expose current hard-coded `max_peak_len=5000`:

```text
--sitega-max-peak-len INT
```

Expected effect:

- shorter peak windows reduce SiteGA scan/evaluation work;
- may change biological sensitivity if peaks are truncated.

This should be documented as a runtime/quality tradeoff.

## Tests

Add unit tests for CLI:

- parser help keeps a single `--jobs` option;
- parser does not include `--bootstrap-jobs`;
- accepts `1`, positive integers, and `-1`;
- rejects `0` and values less than `-1`;
- setup passes resolved positive jobs into `DeNovoPipeline` for bootstrap;
- discovery tool setup passes `threads=1` to threaded discovery tools even when
  `--jobs > 1`;
- MIMOSA comparator setup keeps using its existing internal `n_jobs` path;
- Tomtom comparator setup behavior remains unchanged.

Add unit tests for bootstrap task construction:

- correct number of tasks for a small grid;
- odd/even split indices match current behavior;
- each task has a unique output directory;
- task ordering is stable.

Add unit tests for sequential compatibility:

- `jobs=1` produces the same motif names and statistics as the current
  sequential path using a fake discovery tool.

Add unit tests for process mode:

- use a top-level pickleable fake discovery tool;
- run a tiny grid with `jobs=2`;
- assert the same motif names/statistics as sequential mode;
- assert result order is deterministic.

Add SiteGA seed tests:

- `SitegaDiscoveryTool(seed=123)` passes `seed=123` to `sitega.train`;
- bootstrap task seed derivation produces unique per-task seeds from `--seed`;
- without base seed, SiteGA task seed remains `None` and current behavior is
  preserved.

Avoid real external tools in normal tests. Keep real SiteGA/MEME/STREME runs
behind the existing opt-in external fullrun mechanism.

## Validation commands

Fast validation:

```bash
uv run pytest tests/test_cli_tools.py tests/test_sitega_discovery.py -q
uv run pytest tests/test_pipeline_selection.py tests/test_core_api.py -q
uv run pytest -q
```

Optional SiteGA smoke test after implementation:

```bash
uv run python -c 'import sitega; print(sitega.__file__)'
```

Optional real-data validation should remain opt-in because SiteGA is slow:

```bash
HORDEMOTIFS_RUN_FULLRUN=1 \
HORDEMOTIFS_FULLRUN_LOG_DIR=/tmp/hordemotifs-real-fullrun-logs \
HORDEMOTIFS_STREME_COMMAND=/home/anton/miniconda3/envs/motifhorde/bin/streme \
HORDEMOTIFS_JAVA_COMMAND=/home/anton/miniconda3/bin/java \
/home/anton/miniconda3/envs/motifhorde/bin/python -m pytest \
  tests/test_external_fullrun.py::test_real_data_full_pipeline_smoke_with_verbose_log \
  -q
```

Benchmark recommendation for SiteGA:

```bash
time motifhorde PEAKS.fa PEAKS_gb.fa promoters.fa out-seq \
  --tool sitega \
  -n 10 \
  --lpd 10-40-10 \
  -l 10-16-2 \
  -m auROC \
  -c mimosa \
  --mimosa-metric dice \
  --jobs 1 \
  --seed 42

time motifhorde PEAKS.fa PEAKS_gb.fa promoters.fa out-par \
  --tool sitega \
  -n 10 \
  --lpd 10-40-10 \
  -l 10-16-2 \
  -m auROC \
  -c mimosa \
  --mimosa-metric dice \
  --jobs 6 \
  --seed 42
```

Compare:

- wall time;
- successful completion;
- number of bootstrap motifs;
- no SIGSEGV;
- no duplicate SiteGA seeds in logs when `--seed` is used.

## Risks and mitigations

### Pickling failures

Risk:

- `ProcessPoolExecutor` requires the discovery tool and returned models to be
  pickleable.

Mitigation:

- keep worker functions top-level;
- keep task data as plain dictionaries, strings, integers, and lists;
- test process mode with a top-level fake discovery tool;
- if a specific discovery tool cannot be pickled, add a clear error message and
  keep sequential fallback for that tool.

### Memory pressure

Risk:

- `--jobs N` may launch `N` external tools or JVMs during bootstrap.

Mitigation:

- discovery tool internal threads are set to `1`;
- document conservative `--jobs` values for Java tools;
- keep final discovery sequential in the first patch.

### Duplicate seeds

Risk:

- parallel SiteGA jobs without explicit seeds may share `time(NULL)` values.

Mitigation:

- derive per-task SiteGA seeds from `--seed`;
- document that reproducible parallel SiteGA runs should use `--seed`.

### Noisy child process output

Risk:

- external tools or SiteGA `printf` output may interleave on stdout.

Mitigation:

- keep per-task output directories and log files;
- do not solve stdout capture in the first patch unless interleaving becomes a
  practical problem.

### Behavior changes for discovery tools

Risk:

- tools that previously received `threads=args.jobs` now receive one internal
  thread for discovery.

Mitigation:

- this is intentional to keep `--jobs` from becoming a multiplier;
- document the policy in README;
- benchmark representative tools after the patch.

## Implementation order

1. Keep the CLI surface unchanged: no new jobs option.
2. Update `--jobs` help text to explain bootstrap process parallelism and
   comparator parallelism.
3. Resolve `args.jobs` once for bootstrap worker count with `_resolve_jobs`.
4. Configure discovery tools with one internal thread where supported.
5. Keep comparator setup on its existing path; do not add wrapper-level
   parallelism for MIMOSA or Tomtom.
6. Add `jobs` and `seed` plumbing through CLI setup into `DeNovoPipeline`.
7. Add `jobs` and `seed` to `Bootstrapper`.
8. Extract current odd/even split logic into a small pure helper.
9. Add task construction that writes per-task FASTA files under one temporary
   root.
10. Add top-level discovery worker and result merging.
11. Add sequential execution path using the same task/result code.
12. Add process execution path for `jobs > 1`.
13. Add SiteGA seed support and per-task seed derivation.
14. Add tests for CLI, tool thread policy, comparator jobs, task construction,
    sequential compatibility, process mode, and SiteGA seed passing.
15. Update README with the new `--jobs` semantics and recommended SiteGA
    command.
16. Run fast test suite.
17. Optionally run real-data SiteGA benchmark.

## Acceptance criteria

- No new public jobs parameter is added.
- Existing `--jobs` validation remains simple and consistent.
- `--jobs 1` follows the same discovery/evaluation order as before.
- `--jobs N` runs independent bootstrap discovery tasks in separate processes.
- Discovery tools do not use internal multithreading during bootstrap.
- MIMOSA comparison uses its existing internal `n_jobs` support with the same
  resolved jobs count.
- Tomtom behavior is unchanged except for any existing `n_jobs` value it already
  receives.
- SiteGA bootstrap tasks do not share process memory.
- SiteGA parallel bootstrap with `--seed` uses deterministic unique seeds per
  task.
- Results are returned in deterministic task order.
- Normal test suite passes.
- No real external fullrun is added to default test execution.
