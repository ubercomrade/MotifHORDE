# SiteGA Rewrite and Optimization Plan

This document proposes how to speed up the SiteGA implementation. The key
recommendation is not to keep patching `src/sitega/andy05cell.cpp` in place.
The current file is old C-style C++ with global state, fixed-size arrays,
manual memory management, mixed I/O and computation, and several correctness
risks. Substantial optimization inside that file will make it harder to verify
and maintain.

The safer approach is to keep the current implementation as a reference and
write a new, cleanly designed implementation next to it. The new implementation
must be compared against the current one for both output quality and runtime
before it replaces the existing backend.

## Goals

The rewrite should optimize for:

- correctness first;
- reproducible results with a fixed seed;
- clear data flow and isolated side effects;
- simpler structures instead of global mutable state;
- faster single-threaded execution;
- measurable parity with the current implementation;
- maintainability rather than local micro-optimizations.

The rewrite should not use OpenMP or other parallel execution as the primary
speedup path.

## Why Not Patch the Existing File

`andy05cell.cpp` mixes many responsibilities in one file:

- FASTA parsing;
- sequence encoding;
- background statistics;
- position weighting;
- population initialization;
- mutation and recombination;
- fitness evaluation;
- matrix inversion;
- logging;
- output file generation;
- Python-facing training behavior.

The hot path depends on global arrays and mutable `town` objects. Most candidate
changes copy full per-sequence placement arrays and then recompute fitness from
scratch. This makes local optimization risky: a small change can alter the
search trajectory, change RNG usage, or hide a bug behind faster execution.

The current code should remain available as:

- a reference implementation;
- a regression oracle for fixed seeds and small datasets;
- a fallback backend while the new implementation is validated.

## Proposed Architecture

Create a new SiteGA backend with small, explicit modules:

- input parsing and validation;
- sequence encoding;
- background statistics;
- k-mer position scoring;
- candidate representation;
- fitness evaluation;
- mutation operators;
- recombination operators;
- search loop;
- output formatting.

Keep the public Python API stable. The Python wrapper can later select the new
backend by default after parity and speed are demonstrated.

Suggested C++ structures:

- `EncodedSequences`: owns forward and reverse dinucleotide arrays, lengths,
  and optional prefix-count tables.
- `BackgroundStats`: owns `dav`, `dcv`, `pexp`, and k-mer log-ratio tables.
- `PositionWeights`: owns selected candidate positions and cumulative weights
  per sequence.
- `Feature`: replacement for `uno`, storing `start`, `end`, and dinucleotide
  code.
- `Candidate`: replacement for `town`, storing features, sequence placements,
  cached score components, and a fingerprint.
- `FitnessWorkspace`: reusable temporary buffers for covariance, feature
  vectors, and linear solves.
- `SearchConfig`: simple configuration structure derived from `TrainParams`.

Use plain structs and vectors. Avoid speculative class hierarchies and service
objects.

The current extension is built through `setuptools` and `Pybind11Extension` with
C++17. Do not add CMake, C++20-only features, or a second build system unless a
separate tooling decision justifies that cost. Keep the rewrite compatible with
the existing package build first; build-system migration is outside the critical
path for making SiteGA correct and faster.

The C-facing `sitega_train()` boundary should remain a small status-returning
API. Internal C++ code may use exceptions for unrecoverable validation or I/O
errors, but they must be caught before crossing the C ABI and translated to an
explicit non-zero status and Python-facing failure.

If new options are needed, such as backend selection or `output_count`, add them
as optional fields with zero/default values that preserve current behavior.
Avoid changing the public Python return tuple until there is a clear need.

## Build and Tooling Requirements

Before relying on the new backend, add a debug/test build path for the extension:

- normal release builds should keep using optimized `-O3` settings;
- development builds should enable at least `-Wall`, `-Wextra`, and
  `-Wpedantic` where supported;
- sanitizer builds should be available for AddressSanitizer and
  UndefinedBehaviorSanitizer;
- sanitizer failures and compiler warnings in changed C++ code should be fixed,
  not documented as acceptable;
- OpenMP should remain optional and should not be required for correctness or
  the main speedup claim.

These checks should be wired into the existing setuptools/pybind11 workflow
unless the project later chooses a broader build-system migration.

## Correctness Fixes to Carry Into the Rewrite

The new implementation should intentionally decide how to handle these current
issues:

- In final scoring, `EvalMahFITTrain()` omits `f1[k] /= nseq` before covariance
  and `df` computation. The nearby `df[k] /= nseq` divides a zero-initialized
  array and has no useful effect. This should be fixed and covered by a test.
- The `infc` parameter computes an information-content value, but that value
  is not used in the current final fitness. The rewrite must either preserve
  current behavior explicitly or implement the intended behavior behind a
  controlled option.
- Weighted position selection is weakened by `dw = maxw` where the code likely
  intended `maxw = dw`. The rewrite should preserve this only in
  `legacy_compat`; corrected mode should implement the intended weighting and
  compare the effect.
- Position sampling currently uses forward-coordinate high-scoring windows even
  when `ori == 1`; only the `fpr` lookup mirrors the reverse-strand position
  back to forward coordinates. The rewrite should either preserve this in
  `legacy_compat` or intentionally switch to strand-consistent sampling.
- Recombination duplicate checks should use the actual shuffled parent indices
  from `pair_d[pair_take[k]]`, not stale pair-list indices from `pair_d[k]`.
- `Reco2_Two_dinucleotides_full()` appears to update `odg[]` with group-relative
  indices instead of offsets for the selected dinucleotide groups. This should
  be repaired or intentionally excluded from compatibility-sensitive paths.
- Validate `motif_len`, `max_lpd`, `size`, `olig_bg`, sequence lengths, and
  population size before allocation. Validation should include current fixed
  limits such as `motif_len <= MOTLEN`, `max_lpd <= LPDLEN`,
  `size <= POPSIZE`, `size <= 16 * (motif_len - 1)`, `olig_bg <= motif_len`,
  and enough candidate motif windows to select a valid top-third threshold.
- Track invalid dinucleotide intervals explicitly. Current training fitness
  returns zero if any selected LPD interval contains `-1`; the rewrite should
  preserve or intentionally replace that behavior with tests.
- Replace `exit(1)` inside library code with explicit error returns or
  exceptions caught before the C ABI boundary and translated into Python-visible
  failures.

## Main Speedup Opportunities

Do the speedups in stages. First make a compatible full-fitness implementation
that can be compared against fixed candidates from the old backend. Add
incremental fitness and rollback only after that full recomputation path is
tested, because those optimizations make cache invalidation and search-trajectory
bugs much harder to isolate.

### 1. Incremental Fitness for Placement Mutations

The hottest mutation, `MutRegShiftHoxaW`, changes one sequence placement:

- one sequence index;
- its motif start;
- possibly its strand.

The current code recomputes the full fitness over all sequences after this
change. The rewrite should cache each candidate's contribution sums:

- feature mean sums;
- covariance second-moment sums;
- optional positional information-content counts;
- k-mer score sum.

For a placement mutation, update the cached state by subtracting the old
sequence contribution and adding the new one. This changes the cost from roughly
O(`nseq * feature_count * lpd_length`) to O(`feature_count * lpd_length`) plus
the matrix solve.

This is likely the largest single-threaded speedup.

This should not be the first implementation target. It depends on a correct
candidate state model, a validated full recomputation path, and tests that prove
the cached state matches a fresh recomputation after accepted and rejected
mutations.

### 2. Prefix Counts for LPD Feature Values

Feature values are currently computed by scanning every dinucleotide in each LPD
interval. Precompute per-sequence, per-strand prefix counts for the 16
dinucleotide codes:

```text
count(code, start, end) =
    prefix[code][end + 1] - prefix[code][start]
```

Then each LPD feature value is O(1). This speeds up:

- full fitness recomputation;
- incremental placement updates;
- final location scanning;
- output model scoring.

Memory cost is predictable: `2 * nseq * (len - 1) * 16` counters. If this is too
large for long peaks, use compact integer types or build prefix tables per
sequence.

### 3. Avoid Full Matrix Inversion

Fitness needs:

```text
mah = df^T * covariance^-1 * df
```

It does not need the full inverse matrix. Solve:

```text
covariance * x = df
mah = df^T * x
```

Use a small dense linear solver:

- Cholesky or LDLT when the covariance matrix is positive definite enough;
- a regularized fallback when the matrix is near singular.

This reduces work, improves numerical behavior, and avoids mutating a global
matrix as the API.

### 4. Candidate Fingerprints for Duplicate Checks

The current duplicate check compares every candidate against the population by
scanning all placements and all features. Add stable fingerprints:

- one hash for `features`;
- one hash for `pos/ori`;
- one combined candidate hash.

Use a hash set for population membership. Only perform full equality checks
when hashes collide. This removes repeated O(`population_size * nseq`) scans
from the mutation and recombination hot paths.

Fingerprints are an acceleration structure, not the source of truth. Keep a full
equality fallback for collisions and test that feature ordering, placement
ordering, and strand changes all affect the combined candidate identity.

### 5. In-Place Mutation With Rollback

The current mutation loop copies the entire candidate before each attempt. The
new implementation should mutate in place and keep a small rollback record:

- old feature values for feature mutations;
- old position and strand for placement mutations;
- old cached score components when using incremental fitness.

If the candidate is rejected, restore from the rollback record. If accepted,
commit the cached state and fingerprint.

This avoids repeated allocation and full-array copies.

Rollback should be introduced after the pure mutation operators and duplicate
checks are covered by invariant tests. Until then, full-copy mutation is slower
but easier to debug.

### 6. Limit Final Output Work

The current implementation writes `*_matN` and `*_locN` for every population
member. The Python wrapper usually consumes only the first few motifs.

The new backend should support:

- `output_count`;
- default output count matching the requested motif count;
- optional legacy mode that writes all population members.

This can remove a large amount of final scanning work.

### 7. Faster Position Weight Preprocessing

Current preprocessing sorts all window scores per sequence to find a top-third
threshold. The rewrite can use `nth_element` for threshold selection, then build
the selected-position list in linear time.

This is not the main bottleneck, but it is simple and preserves clarity.

### 8. Remove Unused Work From the Hot Path

If preserving current behavior, do not compute values that do not affect
selection:

- skip information-content computation when `infc` is inactive or when current
  compatibility mode ignores `inf`;
- avoid log/debug formatting unless verbose logging requires it;
- avoid building temporary arrays that are not consumed.

These changes should be controlled by tests because they can alter floating
point operation order or RNG consumption if done carelessly.

## Validation Strategy

The rewrite should be validated in stages.

### Stage 1: Deterministic Component Tests

Add unit tests for pure functions:

- sequence encoding;
- reverse complement encoding;
- background means and covariance;
- k-mer log-ratio computation;
- position candidate selection;
- LPD feature value computation;
- covariance assembly;
- Mahalanobis score for a fixed small example;
- linear-solver parity against the current matrix-inversion result on small
  non-singular examples;
- mutation invariants;
- recombination invariants.

Add a small debug-only way to construct and score a fixed candidate without
running the full genetic search. This gives a stable comparison target for the
new full-fitness implementation before RNG-driven search behavior is involved.

### Stage 2: Golden Tests Against Current Backend

Create small FASTA fixtures and run both implementations with the same seed.
Compare:

- full fitness for fixed candidate states;
- best fitness;
- generated matrix files;
- location files;
- candidate feature invariants;
- ranking stability where exact equality is reasonable.

Exact parity may not hold after fixing known bugs. In that case, keep two
comparison modes:

- `legacy_compat`: reproduces current component behavior as closely as
  practical for validation;
- `corrected`: applies intentional bug fixes and compares biological/output
  quality instead of exact byte-level parity.

Do not make byte-for-byte reproduction of the whole stochastic search a hard
requirement after bug fixes. It is enough to prove component parity where
intended, deterministic behavior for the new backend, and stable output quality
on fixed benchmark inputs.

### Stage 3: Runtime Benchmarks

Benchmark both implementations on:

- tiny test data for CI;
- small real-data subset;
- 400-record real-data subset used by external smoke tests;
- at least one larger dataset if available.

Measure:

- total runtime;
- number of fitness evaluations;
- time spent in fitness evaluation;
- time spent in duplicate checks;
- final output writing time;
- best fitness and output motif count.

Benchmarks should report fixed seed, input sizes, `motif_len`, `size`,
`pop_size`, and output count.

### Stage 4: Integration Rollout

Add the new backend without removing the old one:

- keep old `sitega.train` behavior available;
- expose an opt-in backend selector;
- expose `output_count` only as an optional setting that defaults to legacy full
  population output until the wrapper is ready to request fewer motifs;
- run both backends in CI on small deterministic data;
- switch the default only after speed and result quality are documented.

## Compatibility Policy

The rewrite should distinguish between compatibility and correctness.

Preserve:

- public Python return shape;
- output file format where practical;
- command-line parameters;
- deterministic seed behavior for the new backend;
- C ABI status-return semantics at the `sitega_train()` boundary.

Do not preserve:

- global mutable state;
- accidental bugs unless explicitly required for `legacy_compat`;
- full-population output when the caller asks for fewer motifs;
- exact floating point artifacts from full matrix inversion;
- hidden calls to `exit(1)` from library code.

## Suggested Implementation Order

1. Add a debug/sanitizer build path and stricter warnings for the SiteGA
   extension.
2. Add profiling counters to the current implementation to establish baseline
   runtime and hot-path call counts.
3. Add a fixed-candidate scoring harness for comparing old and new fitness
   implementations without running the full genetic search.
4. Write new sequence encoding and prefix-count code with unit tests.
5. Implement background statistics and position weighting.
6. Implement candidate representation without cached incremental state.
7. Implement full fitness recomputation using prefix counts and linear solve.
8. Compare full fitness on fixed candidates against the old implementation.
9. Implement duplicate checks with fingerprints and full-equality fallback.
10. Implement mutations and recombinations with invariant tests.
11. Add incremental fitness for placement mutations and verify it against fresh
    full recomputation after each mutation type.
12. Add in-place mutation rollback.
13. Implement output generation with configurable output count.
14. Run golden tests and runtime benchmarks against the old backend.
15. Add an opt-in Python/backend switch.
16. Make the new backend default only after validation.

## Expected Impact

The most important expected speedups are:

- large reduction from incremental placement mutations;
- constant-time LPD feature scoring from prefix counts;
- faster and more stable Mahalanobis computation by solving a linear system;
- removal of repeated full population scans through candidate fingerprints;
- less final output work when only top motifs are needed.

Together these changes should be capable of multi-fold single-threaded speedups
without changing the external API or relying on OpenMP.
