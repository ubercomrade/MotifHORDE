# Plan: directed SiteGA feature selection

## Scope

This plan covers only internal feature selection in
`src/sitega/andy05cell_update.cpp`.

Model/search hyperparameters remain external:

- `motif_len`
- `max_lpd`
- `size` / `feature_count`
- `olig_bg`
- `pop_size`
- `generations`
- `mutation_attempts`
- `stale_generations`

The target is to improve how the algorithm chooses model features
`Feature{start, end, code}` and how motif placements are adjusted after those
features change.

The first implementation should keep the existing random mutations and
recombination. The new logic should be added as one extra directed search path,
so behavior remains easy to compare and roll back.

## Current Weak Point

The current implementation has a directed fitness function, but feature
proposals are mostly random:

- `make_random_candidate` samples feature `code`, `start`, and `end` randomly.
- `mutate_feature_code` randomly replaces only the dinucleotide code.
- `mutate_feature_interval` randomly replaces only the interval.
- A feature mutation is evaluated against the candidate's current
  `positions`/`orientations`.

This means a good feature can be rejected when the current placements are poor.
It also means many mutation attempts are spent on features that foreground and
background data already suggest are weak.

## Intended Direction

Add two cooperating mechanisms:

1. Precompute a ranked pool of all possible features from foreground/background
   statistics.
2. Add a directed feature replacement mutation that immediately performs a
   short placement refinement before accepting or rejecting the new feature set.

The combined flow should be:

```text
encoded foreground + background stats
    -> ranked feature pool
    -> replace one candidate feature from the pool
    -> repair/refine placements for the changed feature set
    -> recompute final fitness
    -> accept only if fit improves and candidate is not duplicated
```

## New Data Structures

Add a small score wrapper near `Feature`.

```cpp
struct FeaturePoolEntry {
    Feature feature;
    double score = 0.0;
    double signed_effect = 0.0;
    double fg_mean = 0.0;
    double bg_mean = 0.0;
    double variance = 0.0;
};
```

Keep this as a simple struct. Do not introduce a class, registry, or separate
service object.

Recommended internal constants:

```cpp
constexpr double kFeatureScoreRidge = 1e-9;
constexpr int kDirectedFeatureAttempts = 64;
constexpr int kFeaturePoolTopMultiplier = 12;
constexpr int kPlacementRefineSequenceLimit = 32;
constexpr int kPlacementRefineSamples = 12;
```

These are search-budget constants, not public model hyperparameters. Keep them
local to this backend unless later evidence shows they need to be configurable.

## Feature Pool Construction

Add `build_feature_pool` after `build_background_stats` or near the other
feature helper functions.

Suggested signature:

```cpp
std::vector<FeaturePoolEntry> build_feature_pool(
    const EncodedSequences& sequences,
    const BackgroundStats& background,
    const SearchConfig& config
);
```

Enumeration:

```text
for code in 0..15
for start in 0..motif_len - 2
for end in start..min(start + max_lpd - 1, motif_len - 2)
```

For each feature, estimate foreground statistics over plausible foreground
windows:

- Use each sequence's existing `candidate_positions`.
- Evaluate both orientations.
- Use the same placement coordinate convention as the current candidate code:
  for reverse orientation, `position` is a coordinate on `sequence.reverse`.
- Skip invalid intervals.
- Normalize per sequence before averaging globally, so longer sequences do not
  dominate the pool score.

A practical foreground statistic:

```text
sequence_mean = average feature value across valid candidate placements
fg_mean = average sequence_mean across sequences with at least one valid placement
fg_var = variance of sequence_mean
```

Use background statistics already computed by `build_background_stats`:

```text
span = end - start
bg_mean = background.mean[span][code]
bg_var = background.covariance[span][code]
```

Score:

```text
diff = fg_mean - bg_mean
variance = fg_var + bg_var + kFeatureScoreRidge
score = diff * diff / variance
signed_effect = diff / variance
```

Use squared score for ranking so both enriched and depleted features can be
selected. Keep `signed_effect` for optional placement scoring or diagnostics.

Sorting:

- Sort descending by `score`.
- Drop entries with non-finite score.
- Keep the full sorted pool initially. Limit top-N sampling at use sites rather
  than truncating the pool during construction.

## Shared Feature Value Helper

Add a small pure helper so pool scoring and placement refinement do not need to
construct a temporary `Candidate`.

Suggested helper:

```cpp
bool feature_value_for_placement(
    const Feature& feature,
    const EncodedSequence& sequence,
    int position,
    unsigned char orientation,
    double& value
);
```

Behavior:

- Select the prefix table by `orientation`.
- Compute `start = position + feature.start`.
- Compute `end = position + feature.end`.
- Return `false` if the interval contains invalid bases.
- Otherwise set `value` to the interval frequency of `feature.code` and return
  `true`.

Then `feature_values_for_placement` can reuse this helper. This removes a small
piece of duplication and keeps invalid-interval handling consistent.

## Directed Feature Replacement

Add a new mutation path instead of replacing the existing random mutations.

Suggested signature:

```cpp
bool try_directed_feature_replacement(
    Candidate& candidate,
    const std::vector<Candidate>& population,
    int population_index,
    const std::vector<FeaturePoolEntry>& feature_pool,
    const EncodedSequences& sequences,
    const BackgroundStats& background,
    const SearchConfig& config,
    Rng& rng
);
```

Mutation behavior:

1. Save old `features`, `positions`, `orientations`, and `fingerprint`.
2. Pick one existing feature index at random.
3. Pick a proposed replacement from the top part of the feature pool.
4. Require `can_use_feature(candidate.features, proposed, feature_index)`.
5. Require the replacement to be different from the old feature.
6. Sort features and mark stats invalid.
7. Run short placement repair/refinement.
8. Recompute candidate stats and final fitness.
9. Accept only if fitness improved and the candidate is not a duplicate.
10. Otherwise restore all saved feature and placement state.

Top-pool limit:

```text
pool_limit = min(
    feature_pool.size(),
    max(config.feature_count * kFeaturePoolTopMultiplier, config.feature_count)
)
```

This keeps proposals directed without making every candidate identical.

If no valid replacement is found after `kDirectedFeatureAttempts`, return
`false`.

## Placement Repair And Refinement

Add placement refinement specifically for the directed feature mutation. Do not
run it after every existing random mutation in the first version.

Suggested signature:

```cpp
void refine_placements_after_feature_change(
    Candidate& candidate,
    const EncodedSequences& sequences,
    const BackgroundStats& background,
    const SearchConfig& config,
    Rng& rng
);
```

The function should be bounded. Avoid full scans across all windows for every
mutation attempt.

Recommended strategy:

1. Recompute stats for the candidate after feature replacement.
2. If current placements contain invalid intervals, first try to repair them.
3. If the candidate is valid and has usable weights, refine by model score.
4. Recompute final stats once after changing placements.

### Repair Mode

Repair mode is needed because a new feature can make the current placements
invalid. If the code computes model weights before repair, weights may be all
zero and the mutation will be rejected too early.

For selected sequences:

- Always include sequences whose current placement has invalid intervals.
- Also include a small random subset of other sequences up to
  `kPlacementRefineSequenceLimit`.
- Build a small candidate placement set:
  - current position;
  - `kPlacementRefineSamples` calls to `sample_weighted_position`;
  - both orientations for each sampled position.
- Skip placements where any candidate feature interval is invalid.
- In repair mode, choose the valid placement with the best
  `placement_kmer_weight`.

If no valid placement is found for a sequence, leave it unchanged. The final
fitness will remain zero if invalid intervals remain, and the mutation will be
rejected.

### Model-Score Refinement Mode

After repair, if the candidate is valid:

1. Compute temporary model weights for the changed feature set.
2. For selected sequences, evaluate the sampled positions and both
   orientations.
3. Skip invalid placements.
4. Choose the placement with the highest model score.

Use a local placement scoring helper that accepts precomputed feature values.
Do not call `recompute_candidate_stats` for every trial window.

The refinement can use the existing output scoring idea:

```text
score = sum(weights[i] * feature_value[i])
normalized_score = (score - weights.minimum) / weights.range
```

The normalized score is only used for comparing placements inside one
candidate. The final accept/reject decision must still use `candidate.fit`.

## Integration Points

### `train_impl`

After:

```cpp
auto background = build_background_stats(background_sequences, config.max_lpd);
```

add:

```cpp
auto feature_pool = build_feature_pool(encoded, background, config);
```

Log basic diagnostics:

```text
Feature pool size=<n> best_score=<score>
```

Then pass `feature_pool` into `initialize_population` only if directed
initialization is implemented. Always pass it into `run_search`.

### `run_search`

Change signature to accept:

```cpp
const std::vector<FeaturePoolEntry>& feature_pool
```

Change mutation selection from three types to four types:

```text
0 -> random feature code mutation
1 -> random feature interval mutation
2 -> placement mutation
3 -> directed feature replacement + placement refinement
```

Keep old mutation paths unchanged.

Track accepted directed mutations separately in the log if useful:

```text
Gen <n> Fit <fit> Mut <accepted_random> Dir <accepted_directed> Rec <accepted_recombinations> Delta <delta>
```

If the log format should stay stable, include directed mutations in `Mut` and
only add extra verbose logging when `config.verbose` is true.

### `initialize_population`

First implementation can leave initialization unchanged. This isolates the
effect of the new mutation path.

Second implementation can bias feature initialization from the feature pool:

- choose positions/orientations as today;
- fill features from top-pool samples using `can_use_feature`;
- fall back to current random feature generation if the pool cannot fill the
  candidate.

Do not remove the fallback. It preserves diversity and protects edge cases where
the ranked pool is too narrow.

## Acceptance And Rollback Details

Directed mutation must roll back all state affected by replacement/refinement:

- `features`
- `positions`
- `orientations`
- `fingerprint`
- statistics, by recomputing after restore

Do not try to preserve old `feature_sum`, `second_moment`, `kmer_sum`, `mah`,
`fpr`, and `fit` manually. Recompute after restore. This is simpler and less
error-prone.

Acceptance condition should match existing mutation policy:

```cpp
candidate.fit > old_fit + kScoreEpsilon &&
!duplicate_candidate(candidate, population, population_index)
```

If rejected, restore and recompute:

```cpp
candidate.features = old_features;
candidate.positions = old_positions;
candidate.orientations = old_orientations;
candidate.fingerprint = old_fingerprint;
candidate.stats_valid = false;
recompute_candidate_stats(candidate, sequences, background, config);
```

## Runtime Controls

The main runtime risk is placement refinement. Keep it bounded from the first
implementation.

Recommended limits:

- Use at most `kPlacementRefineSequenceLimit` sequences per directed mutation,
  plus any sequences that must be repaired because their current placement is
  invalid.
- Use at most `kPlacementRefineSamples` sampled positions per selected sequence.
- Evaluate both orientations for each sampled position.
- Recompute full candidate stats once before refinement and once after
  refinement, not per trial placement.

If runtime is still high, reduce directed mutation frequency before reducing
population size or generations. For example, use four mutation types but let
directed replacement happen only on every second attempt or only for the better
half of the population.

## Correctness Risks

### Invalid intervals

New features can introduce invalid intervals for current placements. The
placement refinement must try repair before relying on model weights. If invalid
intervals remain, `evaluate_from_stats` will correctly return zero fitness.

### Reduced diversity

If all candidates draw from the same top-ranked features, the population can
collapse. Mitigations:

- keep old random mutations;
- sample from top-N rather than always picking rank 0;
- keep random initialization in the first implementation;
- keep recombination unchanged.

### Local overfitting to candidate windows

The pool is built from `candidate_positions`, which are already biased by
k-mer enrichment. This is intentional, but it can overemphasize the initial
k-mer filter. Mitigation:

- include both orientations;
- normalize per sequence;
- keep placement mutation and recombination active;
- compare results against the current backend on real and synthetic tests.

### More duplicated scoring code

Avoid adding another separate implementation of feature value extraction. Add
`feature_value_for_placement` and reuse it from pool scoring, refinement, and
the existing per-placement feature evaluation where practical.

## Verification Plan

Fast checks:

```bash
uv run ruff check .
uv run pytest
```

Targeted behavioral checks to add or perform:

- Same seed should produce deterministic output.
- Feature pool should be non-empty for valid configs.
- Feature pool should be sorted by descending finite score.
- Directed mutation should restore the candidate exactly when rejected.
- Directed mutation should not introduce duplicate candidates.
- Candidate features should still satisfy `valid_feature_set`.
- Candidates with unresolved invalid intervals should get zero fitness and be
  rejected.

Integration comparison:

- Run the existing fast SiteGA tests.
- Compare logs before and after the change for:
  - accepted directed mutation count;
  - final best fitness;
  - runtime.

External real-data full-run tests remain opt-in. Do not run them for the first
implementation unless specifically validating real integrations.

## Suggested Implementation Order

1. Add `FeaturePoolEntry` and `feature_value_for_placement`.
2. Implement `build_feature_pool` and log pool size/best score.
3. Add tests or lightweight assertions for pool construction behavior.
4. Add `try_directed_feature_replacement` without placement refinement; verify
   rollback and duplicate handling.
5. Add bounded repair/refinement of placements inside the directed mutation.
6. Wire directed mutation into `run_search` as a fourth mutation type.
7. Run fast tests and inspect logs for acceptance rate and runtime.
8. Only after that, consider feature-pool-biased initialization.

This order keeps the diff reviewable and makes it possible to isolate whether
quality changes come from directed feature proposals or from placement
refinement.
