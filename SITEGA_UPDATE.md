# SiteGA Island Model Update Plan

## Goal

Implement an island model for the SiteGA search in
`src/sitega/andy05cell_update.cpp`.

Target semantics:

- `num_motifs` becomes the number of islands.
- Each island runs its own population search.
- `pop_size` defines the number of individuals per island, not the total number
  of individuals across all islands.
- Each island returns exactly one best candidate.
- Output files are written from the island winners, sorted by fitness.
- The wrapper should pass `num_motifs` as the desired number of islands, so
  `num_motifs == num_islands`.

The implementation should keep the current code simple. Do not introduce a
class hierarchy, scheduler abstraction, generic GA framework, or speculative
parallel execution layer.

## Current Behavior To Replace

Current flow in `train_impl`:

1. Build shared encoded foreground sequences.
2. Build shared background statistics.
3. Build a shared feature pool.
4. Create one population with `config.pop_size` candidates.
5. Run `run_search` once on that population.
6. Sort the final population.
7. Write the top `config.num_motifs` candidates.

This means `num_motifs` currently controls how many top candidates are written
from one mixed population. It does not create independent search niches.

## New Behavior

New flow in `train_impl`:

1. Build shared encoded foreground sequences.
2. Build shared background statistics.
3. Build a shared feature pool.
4. Create `config.num_motifs` independent islands.
5. For each island:
   - initialize a fresh population of `config.pop_size` candidates;
   - run the existing search on that island only;
   - sort the island population;
   - keep only `population.front()` as the island winner.
6. Sort all island winners by `fit`.
7. Write one motif per island winner.

No migration is included in the first implementation. Migration is useful later
only if independent islands produce too many weak local optima.

## Configuration Changes

Keep the public ABI stable:

- Do not add fields to `TrainParams`.
- Do not add a new Python API argument.
- Reinterpret `num_motifs` inside the backend as island count.

Update comments and logs so the new meaning is explicit:

- In `src/sitega/sitega_train.h`, change the `num_motifs` comment from
  "number of best motifs to write" to "number of independent islands / motifs
  to write".
- In `src/sitega/bindings.cpp`, update comments/docstrings that describe
  `num_motifs`.
- In `train_impl`, log both:
  - `pop_size_per_island=<config.pop_size>`;
  - `islands=<config.num_motifs>`.

Keep the existing validation cap:

- `pop_size` remains limited by `kMaxPopulation`.
- `num_motifs` remains clamped to `1..pop_size` initially if preserving current
  behavior is important.

Consider a small follow-up cleanup after the island change:

- If `num_motifs` is now independent from `pop_size`, the cap by `pop_size` is
  no longer conceptually required.
- A clearer cap would be a new `kMaxIslands`, for example `100`.
- Do not combine this semantic cleanup with the first implementation unless
  tests make the current cap awkward.

## C++ Implementation Steps

### 1. Add a focused island helper

Add a small function near `run_search` or near `train_impl`:

```cpp
std::vector<Candidate> run_island_searches(
    const std::vector<FeaturePoolEntry>& feature_pool,
    const EncodedSequences& encoded,
    const BackgroundStats& background,
    const SearchConfig& config,
    Rng& rng,
    std::ostream& log
)
```

Responsibility:

- run all islands;
- return one winner per island;
- keep all I/O limited to logging;
- avoid writing output files;
- avoid rebuilding shared data.

Pseudo-flow:

```cpp
std::vector<Candidate> winners;
winners.reserve(static_cast<std::size_t>(config.num_motifs));

for (int island = 0; island < config.num_motifs; ++island) {
    log << "Island " << (island + 1) << "/" << config.num_motifs << '\n';

    auto population = initialize_population(
        encoded,
        background,
        config,
        rng,
        log
    );
    run_search(population, feature_pool, encoded, background, config, rng, log);
    sort_population(population);

    winners.push_back(std::move(population.front()));
    log << "Island " << (island + 1)
        << " winner_fit=" << winners.back().fit << '\n';
}

sort_population(winners);
return winners;
```

Keep this helper simple. It should not own configuration, output paths, or
shared preprocessing.

### 2. Update `train_impl`

Replace:

```cpp
Rng rng(config.seed);
auto population = initialize_population(encoded, background, config, rng, log);
run_search(population, feature_pool, encoded, background, config, rng, log);
sort_population(population);
write_outputs(population, encoded, background, config, result);
```

with:

```cpp
Rng rng(config.seed);
auto winners = run_island_searches(
    feature_pool,
    encoded,
    background,
    config,
    rng,
    log
);
write_outputs(winners, encoded, background, config, result);
```

`write_outputs` can keep accepting `std::vector<Candidate>& population`; rename
the parameter to `candidates` if editing that function anyway.

### 3. Keep RNG deterministic

Use the single existing `Rng rng(config.seed)` and pass it through islands in
sequence.

This preserves deterministic behavior for a fixed seed while ensuring each
island receives a different segment of the RNG stream.

Do not seed every island with `config.seed + island` in the first version. That
looks simple but creates avoidable questions around seed collisions,
reproducibility, and platform-specific integer behavior.

### 4. Keep shared preprocessing outside islands

The following must remain outside the island loop:

- FASTA parsing;
- k-mer log-ratio calculation;
- sequence encoding;
- background statistics;
- feature pool construction.

Only population initialization and evolutionary search should run per island.
This keeps the implementation correct and avoids repeated expensive work.

### 5. Preserve output naming

Continue writing:

- `<base>_mat1`, `<base>_loc1`;
- `<base>_mat2`, `<base>_loc2`;
- ...

The files now correspond to sorted island winners, not the top candidates from
one population.

Keep `TrainResult` unchanged:

- `mat_path` points to `_mat1`;
- `loc_path` points to `_loc1`;
- `best_fit` is the best winner fitness.

## Python Wrapper Changes

Current wrapper call:

```python
num_motifs=max(self.num_motifs, number_of_motifs)
```

This already mostly matches the desired island semantics: request at least as
many islands as motifs requested by the caller.

Recommended wrapper behavior:

```python
num_motifs=max(self.num_motifs, number_of_motifs)
```

Keep this unless the product-level intent is to make `number_of_motifs` the
exact island count. If exact behavior is desired, change it to:

```python
num_motifs=number_of_motifs
```

Given the user's requirement that `num_motifs == num_islands`, the important
backend condition is that the value passed to `sitega.train(..., num_motifs=...)`
is interpreted as island count.

Update comments in `SitegaDiscoveryTool`:

- `pop_size`: individuals per island;
- `num_motifs`: default number of islands / SiteGA motifs to request.

## Test Plan

### Unit-level C++ behavior

If there are existing C++ tests for `sitega_train`, extend them. If not, prefer
Python-level tests first to avoid adding heavy C++ test infrastructure.

The core behavior to verify:

- with `num_motifs = 3`, exactly three `_matN`/`_locN` output pairs are written;
- each output corresponds to one island winner;
- `result.best_fit` equals the best winner fitness;
- fixed seed produces deterministic outputs;
- `pop_size` controls each island population size.

The last point can be checked through logs first:

- log line includes `islands=3`;
- each island logs initialization/search;
- each island reports a winner.

### Python wrapper tests

Update existing tests in `tests/test_sitega_discovery.py`:

- Rename expectations/comments so `num_motifs` means islands.
- Keep the existing assertion that the wrapper passes
  `max(self.num_motifs, number_of_motifs)` if that behavior remains.
- If exact requested motif count is chosen instead, update
  `test_sitega_discovery_requests_at_least_requested_motifs` accordingly.

### Fast integration test

Run the normal fast test suite:

```bash
uv run pytest tests/test_sitega_discovery.py -q
```

If the backend has direct Python tests that compile/import `sitega`, run the
smallest relevant subset as well.

### External full-run test

Do not run the full real-data smoke test by default. It is slow and explicitly
opt-in.

Run it only if the implementation changes are finished and real external
validation is required:

```bash
HORDEMOTIFS_RUN_FULLRUN=1 \
HORDEMOTIFS_FULLRUN_LOG_DIR=/tmp/hordemotifs-real-fullrun-logs \
HORDEMOTIFS_STREME_COMMAND=/home/anton/miniconda3/envs/motifhorde/bin/streme \
HORDEMOTIFS_JAVA_COMMAND=/home/anton/miniconda3/bin/java \
/home/anton/miniconda3/envs/motifhorde/bin/python -m pytest \
  tests/test_external_fullrun.py::test_real_data_full_pipeline_smoke_with_verbose_log \
  -q
```

## Expected Behavioral Changes

Positive changes:

- multiple local optima are more likely to appear in output;
- one strong optimum has less opportunity to dominate all results;
- `num_motifs` now maps directly to independent search attempts;
- `pop_size` becomes easier to reason about as per-island search strength.

Tradeoffs:

- runtime scales roughly linearly with `num_motifs`;
- total candidate evaluations become approximately
  `num_motifs * pop_size * generations * mutation_attempts`;
- memory remains bounded by one island population at a time if winners are kept
  and full island populations are discarded;
- output motifs can still be similar if different islands converge to the same
  basin.

## Optional Follow-up: Diverse Winner Filtering

Island model increases the chance of finding multiple optima, but it does not
guarantee distinct motifs.

After the basic island implementation, consider adding a simple diversity pass
over island winners:

1. Sort winners by `fit`.
2. Select the best winner.
3. Continue scanning winners.
4. Keep a winner only if it is sufficiently different from already selected
   winners.

This requires a candidate distance function. Start with a simple metric:

- feature distance based on `(code, start, end)` differences;
- placement distance based on position/orientation differences per sequence.

Do not add this in the first island-model patch unless duplicate island winners
are a demonstrated problem.

## Risks And Compatibility

- The meaning of `num_motifs` changes from "top-N from one population" to
  "N independent island winners". Update comments and logs to avoid ambiguity.
- Runtime may increase substantially for existing callers that use the default
  `num_motifs = 20`.
- The current default wrapper settings may become expensive:
  `num_motifs=20`, `pop_size=30`, `generations=50`,
  `mutation_attempts=50`.
- Consider lowering wrapper defaults only after measuring runtime. Do not hide
  this change inside the backend patch.
- Since `TrainParams` stays unchanged, ABI compatibility is preserved.

## Definition Of Done

- `num_motifs` is interpreted as island count in the backend.
- `pop_size` is used as population size per island.
- Each island runs independently and returns one winner.
- Output files are written from sorted island winners.
- Shared preprocessing is not repeated per island.
- Logs clearly describe island count and per-island population size.
- Header/binding/wrapper comments no longer describe old top-N semantics.
- Fast tests pass.
- No unnecessary classes, wrappers, global state, or speculative abstractions are
  introduced.
