# Gradual MIMOSA Integration Plan

## Goal

Make `mimosa` the canonical model runtime used by `motifhorde`, while keeping
`motifhorde` focused on discovery orchestration, validation, parameter
selection, external tool execution, and output layout.

The migration must be gradual. The project should keep working after each
phase, and every behavior change must be covered by tests before duplicated
local code is removed.

## Current State

`motifhorde` currently contains local model-runtime code that overlaps with
`mimosa`:

- `src/motifhorde/models.py` defines `GenericModel`, model registry dispatch,
  model readers and writers, scanning, site extraction, PFM reconstruction, and
  threshold calibration.
- `src/motifhorde/batches.py` defines dense `TypedDict`-based sequence,
  score, and strand-profile payloads.
- `src/motifhorde/functions.py` contains shared numerical scoring and
  transformation utilities.
- `src/motifhorde/comparison.py` contains matrix and profile comparison logic
  that is conceptually close to `mimosa`.
- `src/motifhorde/discovery.py`, `src/motifhorde/evaluation.py`, and
  `src/motifhorde/pipeline.py` depend on the local model API.

This makes the project self-contained, but it creates a long-term maintenance
problem: fixes and behavior changes in `mimosa` and `motifhorde` can diverge.

## Target Architecture

The final structure should be:

- `mimosa`: model types, model readers, model writers, scanning, site
  extraction, PFM reconstruction, threshold calibration, and comparison
  algorithms that are not specific to motif discovery orchestration.
- `motifhorde`: de novo discovery wrappers, odd/even validation, performance
  evaluation, parameter-grid handling, final motif selection, output
  persistence, and CLI UX.
- `motifhorde.models`: temporary compatibility module that re-exports or wraps
  the relevant `mimosa` API during the migration.
- `motifhorde.io`: keep only FASTA, pipeline output, and discovery-specific I/O
  that does not belong in the generic model runtime.

The dependency direction must be one-way:

```text
motifhorde -> mimosa
```

`mimosa` must not import from `motifhorde`.

## Guiding Rules

- Do not replace imports blindly. First lock down behavior with tests.
- Prefer moving stable generic model-runtime logic into `mimosa` over keeping a
  fork in `motifhorde`.
- Keep `motifhorde` APIs stable during migration where practical.
- Do not introduce adapter classes unless a plain function is insufficient.
- Keep payload shapes explicit. If adapters are needed, make conversions local
  and easy to delete later.
- Do not mix old and new runtimes in the same execution path without a clear
  compatibility boundary.
- Remove local duplicate code only after tests prove the `mimosa` path is
  behaviorally equivalent.

## Migration Progress

- [x] 2026-06-22: Phase 0 baseline audit recorded for the current local
  runtime. No production behavior was changed.
- [x] 2026-06-22: Phase 1 contract tests added in
  `tests/test_model_contracts.py` for the current `motifhorde.models` API,
  edge cases, pickle/PFM write compatibility, and public comparison wrappers.
- [x] 2026-06-22: `mimosa-tool>=1.3.0,<2` added to project dependencies and
  locked in `uv.lock`.
- [x] 2026-06-22: `src/motifhorde/models.py` converted into a thin
  compatibility facade over `mimosa`, with legacy wrappers only for
  motifhorde-specific defaults and PWM/MEME reader compatibility.
- [ ] Phase 5 and later remain pending. Internal consumers still use the
  compatibility boundary, and local duplicate `batches.py`, `functions.py`,
  `comparison.py`, and compatibility I/O code have not been removed.

## Phase 0 - Baseline Audit

Purpose: define exactly what must stay compatible before any dependency change.

Tasks:

1. Inventory public API currently exported by `motifhorde`.
   - `motifhorde.__init__`
   - `motifhorde.models`
   - documented README examples
   - CLI-visible behavior

2. Inventory local model-runtime features.
   - `GenericModel` fields and pickle behavior
   - `read_model`
   - `write_model`
   - `scan_model`
   - `scan_model_strands`
   - `get_scores`
   - `get_frequencies`
   - `calculate_threshold_table`
   - `get_sites`
   - `get_pfm`
   - `register_model_handler`

3. Map every direct dependency on `motifhorde.models`.
   - `evaluation.py`
   - `discovery.py`
   - `pipeline.py`
   - `comparison.py`
   - `cache.py`
   - tests
   - public imports from `__init__.py`

4. Record differences between local `motifhorde` behavior and current
   `mimosa` behavior.
   - dense `TypedDict` batches vs any `mimosa` sequence payload
   - strand handling: `best`, `+`, `-`, and both-strand output
   - score normalization semantics
   - threshold table format
   - site table column names and dtypes
   - PFM reconstruction pseudocount behavior
   - pickle compatibility
   - supported file extensions such as `.txt` for PWM/MEME input

Deliverable:

- A short compatibility matrix committed as a section in this file or a
  dedicated developer note.

Acceptance criteria:

- No code behavior changes.
- All current tests pass.
- The team knows which local behaviors must be preserved.

Status:

- [x] Completed locally on 2026-06-22.
- [x] Public API and direct internal dependencies are recorded below.
- [x] Compatibility-sensitive runtime behavior is covered by
  `tests/test_model_contracts.py`.

### Compatibility Matrix

| Area | Current `motifhorde` behavior | Migration requirement |
| :--- | :--- | :--- |
| Public API | `motifhorde.__init__` exports `GenericModel`, `get_pfm`, `get_sites`, `read_model`, `scan_model`, and `write_model`. README examples import the fuller `motifhorde.models` API, including `scan_model_strands` and `calculate_threshold_table`. | Preserve these imports through `motifhorde.models` until public docs are intentionally changed. |
| Model container | `GenericModel` is a `motifhorde.models` dataclass with `type_key`, `name`, `representation`, `length`, and `config`. Existing joblib pickles use the `motifhorde.models.GenericModel` class path. | Either load old pickles through a compatibility path or reject them with a clear migration error. |
| Model readers | `pwm` accepts `.meme`, `.txt`, `.pfm`, and compatible `.pkl`; `bamm` accepts `.ihbcp` and resolves a missing suffix; `sitega` accepts `.mat` and `.pkl`; `dimont` and `slim` accept `.xml` and `.pkl`; `scores` accepts FASTA-like numeric profiles. | `mimosa` readers or wrapper code must preserve these accepted inputs before local readers are removed. |
| Dense payloads | Sequences, masked scores, and strand profiles are `TypedDict` payloads with dense arrays, row lengths, masks where needed, and padding values. | Internal consumers must keep receiving the same shapes until they are migrated together. |
| Strand semantics | `scan_model(..., strand="+")` scans the forward strand; `"-"` scans reverse complement; `"best"` returns the elementwise maximum of plus/minus scores. `scan_model_strands` returns a two-strand bundle with shape `(2, n_rows, width)`. | Preserve score values, row lengths, and plus/minus ordering exactly. |
| Threshold table | `calculate_threshold_table` uses valid scores only and returns descending `[score, -log10(tail_probability)]` rows. Empty calibration produces `[[0.0, 0.0]]`. | Preserve lookup semantics and empty input behavior. |
| Site table | `get_sites` returns columns `seq_index`, `start`, `end`, `strand`, `score`, `log_tail`, and `site`. Threshold mode requires `fpr_threshold`; best mode returns one hit per sequence when available. | Preserve column names, ordering, strand labels, and threshold validation. |
| PFM reconstruction | `get_pfm` reconstructs a normalized `4 x length` PFM from selected sites with default pseudocount `0.25`. `top_fraction` keeps at least one site. Empty selections raise `ValueError("No sites found")`. | Preserve normalization, pseudocount behavior, and error behavior. |
| Comparison | `TomtomComparator` and `UniversalMotifComparator` are pipeline-facing wrappers. Heterogeneous motif comparison falls back to PFM reconstruction when `type_key` differs. One-to-many comparison preserves target order. | Delegate only after result parity is covered by contract tests. |
| Direct internal dependencies | `evaluation.py`, `discovery.py`, `pipeline.py`, `comparison.py`, `cache.py`, `io.py`, tests, and `__init__.py` directly import or depend on `motifhorde.models`. | Migrate consumers through one compatibility boundary rather than mixing local and `mimosa` runtimes. |
| Current `mimosa` state | `mimosa-tool 1.3.0` is installed and provides `GenericModel`, dense batches, scanning, both-strand profiles, sites/PFM reconstruction, threshold calibration, readers/writers, and comparison API. | Keep the `motifhorde.models` facade until legacy `.txt`/compact MEME parsing, old defaults, and pickle compatibility are intentionally retired or fully upstreamed. |

## Phase 1 - Contract Tests Before Integration

Purpose: make runtime behavior measurable before changing implementation.

Tasks:

1. Add model API contract tests.
   - `read_model` returns the expected model fields.
   - `scan_model(..., strand="best")` returns stable values and masks.
   - `scan_model(..., strand="+")` and `scan_model(..., strand="-")` are
     shape-compatible.
   - `scan_model_strands` returns plus and minus profiles with identical row
     lengths.
   - `get_sites` returns stable columns: `seq_index`, `start`, `end`, `strand`,
     `score`, `log_tail`, `site`.
   - `get_pfm` returns a normalized 4 x L matrix.

2. Cover each supported model family.
   - `pwm`
   - `bamm`
   - `sitega`
   - `dimont`
   - `slim`
   - `scores`

3. Add edge-case tests.
   - empty sequence batch
   - sequence shorter than motif
   - threshold mode without threshold
   - threshold mode with explicit threshold table
   - top-fraction PFM reconstruction
   - pickle load for existing `motifhorde.models.GenericModel`

4. Add comparison contract tests.
   - same-type motif comparison
   - heterogeneous comparison with PFM reconstruction
   - profile comparison with fixed seed
   - one-to-many comparison result ordering

Deliverable:

- A focused test suite that describes current behavior.

Acceptance criteria:

- Tests fail if scan values, masks, site tables, PFM reconstruction, or
  comparison semantics change unexpectedly.
- Existing full pipeline tests still pass.

Status:

- [x] Completed for the local baseline on 2026-06-22 in
  `tests/test_model_contracts.py`.
- [x] The suite covers exact PWM scan values, masks, strand bundle layout,
  threshold tables, site-table schema and values, PFM reconstruction,
  empty/short sequence batches, threshold validation, supported model-family
  loading, current `GenericModel` pickle loading, PWM write/read roundtrip,
  same-type comparison, heterogeneous PFM-based comparison, fixed-seed profile
  comparison, and one-to-many target ordering.
- [x] The suite now also covers the installed `mimosa` import/API smoke test,
  `motifhorde.models.GenericModel` aliasing to `mimosa.GenericModel`, legacy
  PWM `.txt` input, and old `motifhorde.models.GenericModel` pickle loading
  through the compatibility facade.
- [ ] Real XML fixtures for `dimont` and `slim` are still needed before local
  XML readers can be safely deleted. The current baseline covers their
  supported `.pkl` loading path and shared scanner behavior.

## Phase 2 - Upstream Generic Improvements To MIMOSA

Purpose: make `mimosa` capable of replacing the local runtime without losing
current `motifhorde` behavior.

Tasks:

1. Move or reimplement generic dense-batch functionality in `mimosa`.
   - sequence batch representation
   - masked score batch representation
   - strand profile bundle representation
   - row slicing helpers
   - valid-value flattening helpers

2. Move or reimplement shared numerical kernels in `mimosa`.
   - batch scoring for tensor models
   - both-strand scanning in one backend call where available
   - empirical log-tail table construction
   - score-to-tail lookup
   - PFM/PWM conversion helpers

3. Align `mimosa` model handlers with local behavior.
   - PWM reader supports `.meme`, `.txt`, `.pfm`, and compatible `.pkl`.
   - PWM stores the source PFM under the agreed key.
   - BaMM reader resolves `.ihbcp` consistently.
   - Dimont and Slim XML readers match local tensor shapes.
   - SiteGA reader preserves min/max score bounds.
   - Scores reader returns the same masked profile semantics.

4. Add `scan_model_strands` to `mimosa` if not already present with compatible
   behavior.

5. Add `mimosa` tests equivalent to the new `motifhorde` contract tests.

Deliverable:

- A `mimosa` release or pinned commit that provides the required runtime API.

Acceptance criteria:

- `mimosa` can run the contract tests against the same test fixtures.
- No `mimosa` code imports from `motifhorde`.
- `mimosa` exposes a small, documented public API rather than relying on
  private modules.

Status:

- [x] `mimosa-tool 1.3.0` provides the required dense batch payloads,
  tensor-model scan kernels, both-strand scanning, threshold calibration,
  site/PFM reconstruction, model handler registry, and `scan_model_strands`
  behavior needed by the current `motifhorde.models` facade.
- [x] Installed `mimosa` API was checked by contract tests through the
  `motifhorde.models` compatibility boundary.
- [ ] Upstream parity is not complete: PWM `.txt` files and compact MEME
  headers such as `w=6` are still handled by the local compatibility wrapper.
  Keep this wrapper until `mimosa` accepts those inputs or `motifhorde`
  intentionally drops that legacy behavior.

## Phase 3 - Add MIMOSA As A Dependency

Purpose: introduce the dependency without changing runtime behavior yet.

Tasks:

1. Add `mimosa-tool` to `pyproject.toml`.
   - Prefer a bounded version range, for example `mimosa-tool>=1.2,<2`, once
     the required API is released.
   - If the needed API is not released, pin a specific commit only temporarily.

2. Update lock files and environment files.
   - `uv.lock`
   - `environment.yml`, if the conda environment remains supported

3. Add an import smoke test.
   - Import `mimosa`.
   - Verify the expected public functions exist.
   - Do not yet route production code through `mimosa`.

Deliverable:

- Dependency added with no behavioral migration.

Acceptance criteria:

- `uv sync` succeeds.
- Current test suite still passes.
- No production path depends on both local and `mimosa` runtime at the same
  time.

Status:

- [x] `mimosa-tool>=1.3.0,<2` is present in `pyproject.toml`.
- [x] `uv.lock` includes `mimosa-tool 1.3.0` and its direct dependency
  additions.
- [x] `tests/test_model_contracts.py` includes an import/API smoke test for
  the installed `mimosa` package.

## Phase 4 - Build A Thin Compatibility Boundary

Purpose: switch imports gradually while preserving the current `motifhorde`
public API.

Tasks:

1. Convert `src/motifhorde/models.py` into a compatibility module.
   - Re-export `mimosa.GenericModel`.
   - Re-export `mimosa.StrandMode`.
   - Re-export or wrap `read_model`, `write_model`, `scan_model`,
     `scan_model_strands`, `get_scores`, `get_frequencies`,
     `calculate_threshold_table`, `get_sites`, `get_pfm`, and
     `register_model_handler`.

2. Keep wrappers only where `motifhorde` has stricter or older behavior.
   - old pickle migration
   - legacy config key aliases
   - README-documented return schemas
   - temporary dense-batch conversion if payloads differ

3. Add explicit migration tests for wrapper behavior.
   - `from motifhorde.models import GenericModel` still works.
   - existing pipeline code can consume a `mimosa` model.
   - old pickled `motifhorde.models.GenericModel` files load or fail with a
     clear migration message.

Deliverable:

- `motifhorde.models` becomes a small compatibility layer.

Acceptance criteria:

- Existing public imports still work.
- Contract tests pass through the `motifhorde.models` API.
- Wrapper code is small and contains no duplicated numerical kernels.

Status:

- [x] `src/motifhorde/models.py` now re-exports/delegates to `mimosa`
  `GenericModel`, model registry, reading/writing, scanning, both-strand
  profiles, score/frequency helpers, threshold calibration, sites, and PFM
  reconstruction.
- [x] Legacy wrappers are limited to compatibility behavior:
  `calculate_threshold_table` keeps the old `strand="best"` default,
  `get_sites`/`get_pfm` preserve old calibration and hit-selection defaults,
  and PWM `.meme`/`.txt` loading keeps motifhorde's permissive MEME parser.
- [x] Migration tests cover `GenericModel` aliasing, old
  `motifhorde.models.GenericModel` pickle loading, legacy `.txt` PWM loading,
  and contract behavior after routing through `mimosa`.
- [ ] Internal consumers still need a Phase 5 audit before any remaining local
  duplicate generic code is deleted.

## Phase 5 - Migrate Internal Consumers

Purpose: make internal code use the compatibility boundary consistently before
removing old implementation code.

Tasks:

1. Update `evaluation.py`.
   - Use the compatibility model API.
   - Keep evaluation-specific logic local.
   - Do not duplicate scan or threshold logic.

2. Update `discovery.py`.
   - Discovery tools should produce `mimosa`-compatible models.
   - External tool execution stays in `motifhorde`.
   - Model parsing should use `mimosa.read_model` through the compatibility
     boundary unless parsing is discovery-specific.

3. Update `pipeline.py`.
   - Keep selection and output policy local.
   - Use `get_pfm` and model writing through the compatibility boundary.
   - Preserve output directory layout and statistics JSON schema.

4. Update `cache.py`.
   - Decide whether model fingerprinting belongs in `mimosa` or remains local.
   - Keep cache key stability for existing workflows where practical.

5. Update `comparison.py`.
   - If `mimosa` exposes equivalent comparison functions, delegate to them.
   - Keep `TomtomComparator` and `UniversalMotifComparator` as user-facing
     `motifhorde` wrappers only if the pipeline needs that API.
   - Remove duplicated comparison kernels only after result parity tests pass.

Deliverable:

- Internal code consistently depends on the compatibility model API or direct
  `mimosa` public API, not local duplicate internals.

Acceptance criteria:

- Pipeline tests pass for PWM, BaMM, SiteGA, Dimont, Slim, and scores where
  fixtures exist.
- Full-run external smoke tests pass when enabled.
- Comparison results stay within documented numerical tolerances.

## Phase 6 - Pickle And File Compatibility Migration

Purpose: avoid silently breaking existing saved models.

Tasks:

1. Decide the supported compatibility policy.
   - Option A: support old `motifhorde` pickles through a migration loader.
   - Option B: reject old pickles with a clear error and provide a conversion
     command.
   - Option C: support both for one minor release, then remove old pickle
     support.

2. Prefer stable serialization for future outputs.
   - Keep `joblib` only if class-path compatibility is acceptable.
   - Consider a simple explicit payload format for model metadata plus arrays.
   - Keep MEME/PFM export for downstream PWM-oriented tools.

3. Add tests for old and new saved models.
   - load old local pickles
   - load new `mimosa` pickles
   - write final models from the pipeline
   - re-read written final models

Deliverable:

- A documented compatibility policy for saved model files.

Acceptance criteria:

- Users do not get silent misloads.
- Compatibility failures include actionable error messages.
- Future output format is documented in README.

## Phase 7 - Remove Local Duplicate Runtime Code

Purpose: finish the migration and reduce maintenance burden.

Tasks:

1. Remove local model-runtime code made obsolete by `mimosa`.
   - model handler registry
   - duplicated scan dispatch
   - duplicated model readers and writers
   - duplicated site/PFM reconstruction
   - duplicated threshold calibration

2. Remove or shrink local batch and numerical utility modules.
   - Keep only functions still specific to `motifhorde`.
   - Move generic helpers to `mimosa`.
   - Delete dead imports and stale tests.

3. Update documentation.
   - README supported model section
   - installation dependencies
   - Python API examples
   - saved model compatibility notes

4. Run cleanup checks.
   - `ruff`
   - unit tests
   - full-run smoke tests where external tools are available

Deliverable:

- `motifhorde` no longer contains a fork of the generic model runtime.

Acceptance criteria:

- Diff removes more duplicate code than it adds.
- `motifhorde` remains structurally simpler.
- Public behavior remains documented and tested.

Status:

- [x] Duplicate model registry, scan dispatch, site/PFM reconstruction, and
  threshold-calibration code was removed from `src/motifhorde/models.py`.
- [ ] Duplicate generic helpers still remain in `src/motifhorde/batches.py`,
  `src/motifhorde/functions.py`, `src/motifhorde/comparison.py`, and the
  compatibility MEME reader path. Do not delete them until Phase 5/6 parity
  and external fixture coverage are complete.

## Phase 8 - Stabilize Release Boundary

Purpose: make future divergence less likely.

Tasks:

1. Add dependency policy.
   - State which `mimosa` version range is supported.
   - Pin upper bounds when public API stability is not guaranteed.
   - Update only through tested dependency bumps.

2. Add cross-project contract tests.
   - Run a minimal test set against the supported `mimosa` release.
   - Include model loading, scanning, PFM reconstruction, and comparison.

3. Add release checklist.
   - bump `mimosa` first when model-runtime behavior changes
   - run `mimosa` tests
   - bump `motifhorde`
   - run `motifhorde` contract and pipeline tests
   - update README if public behavior changes

Deliverable:

- A stable maintenance process for both packages.

Acceptance criteria:

- Future model-runtime fixes happen in one place.
- `motifhorde` does not reintroduce local forks of generic `mimosa` code.

## Suggested Implementation Order

1. Add contract tests in `motifhorde` for the current local behavior.
2. Port dense batch and both-strand runtime support to `mimosa`.
3. Release or pin a `mimosa` version with the required API.
4. Add `mimosa-tool` as a dependency without changing runtime behavior.
5. Turn `motifhorde.models` into a compatibility module.
6. Migrate internal consumers module by module.
7. Migrate comparison logic after model scanning parity is proven.
8. Define pickle compatibility and future serialization policy.
9. Remove duplicate local runtime code.
10. Update README and release notes.

## Main Risks

- Pickle incompatibility because class paths differ between
  `motifhorde.models.GenericModel` and `mimosa.models.GenericModel`.
- Numerical drift in scanning or PFM reconstruction due to different batch
  representations or threshold calibration.
- API mismatch around strand handling, especially `scan_model_strands` and
  both-strand calibration.
- Hidden dependency on local dense `TypedDict` payloads in pipeline and
  comparison code.
- External discovery tools may produce outputs that local parsers accept but
  `mimosa` parsers reject.
- Removing local code before contract tests are complete could create subtle
  selection regressions.

## Rollback Plan

Each phase should be merged independently. If a migration phase fails:

1. Revert only that phase.
2. Keep contract tests that exposed the failure.
3. Fix the missing behavior in `mimosa`.
4. Retry the `motifhorde` integration after the `mimosa` API is corrected.

Do not keep parallel local and `mimosa` runtimes active as a permanent fallback.
That would preserve the original maintenance problem.

## Definition Of Done

The migration is complete when:

- `motifhorde` depends on `mimosa` for generic model runtime behavior.
- `motifhorde` no longer duplicates model scanning, site extraction, PFM
  reconstruction, threshold calibration, or generic comparison kernels.
- `motifhorde.models` is either a small compatibility facade or removed from
  public documentation.
- Old saved model compatibility is either supported or explicitly rejected with
  a documented migration path.
- Contract, comparison, pipeline, and full-run smoke tests pass.
- README describes the new dependency boundary clearly.
