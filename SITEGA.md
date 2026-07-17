# SiteGA Algorithm Notes

This document describes the legacy SiteGA implementation retained in
`src/sitega/__andy05cell_old.cpp`. The default backend is the rewritten
`src/sitega/andy05cell_update.cpp`, whose data structures and search operators
differ from the legacy code described below.

## Main Data Model

The implementation trains a motif model represented as a set of locally
positioned dinucleotide features.

Input sequences are converted to dinucleotide codes:

- `a=0`, `c=1`, `g=2`, `t=3`;
- each adjacent pair is encoded as `4 * first + second`;
- valid dinucleotides therefore have codes `0..15`;
- invalid or ambiguous dinucleotides are stored as `-1`;
- both forward and reverse-complement strands are stored.

Important structures:

- `ss s[16]`: lookup table from dinucleotide code to a text label such as
  `ac`, `gt`, etc.
- `uno`: one LPD feature. It stores `sta`, `end`, and `num`, where `num` is
  the dinucleotide code and `[sta, end]` is the local interval inside the motif.
- `town`: one candidate motif model, also called an individual. It stores:
  `tot[]` LPD features; `deg[16]` counts per dinucleotide code; `odg[]` counts
  by LPD length; per-sequence motif position `pos[]`; per-sequence strand
  `ori[]`; fitness fields `mah`, `fpr`, `inf`, and `fit`; mutation and
  recombination counters.
- `pop[MEGE]`: the population of candidate models.
- `det1`, `det2[2]`: temporary candidate buffers used during initialization,
  mutation, and recombination.
- `uw[POPSIZE][POPSIZE]`: shared work matrix used while computing covariance
  and Mahalanobis score.
- `RNG g_rng`: deterministic xorshift RNG seeded from `TrainParams.seed` or
  `time(NULL)`.

Important derived arrays:

- `seq_real[2][nseq][len-1]`: foreground dinucleotide codes for both strands.
- `dav[reg_max][16]`: background mean frequency for each LPD span index
  `end - sta` and dinucleotide. The actual LPD length is `end - sta + 1`.
- `dcv[reg_max][16]`: background self-covariance for each LPD span index
  `end - sta` and dinucleotide.
- `pexp[16]`: expected dinucleotide frequencies from background.
- `octa_pro1`: per-position foreground/background k-mer log-ratio before motif
  window averaging.
- `octa_prow`: per-position averaged k-mer log-ratio for a full motif window.
- `octa_prowb`: selected high-scoring motif-window positions per sequence.
- `octa_prows`: cumulative weights for selecting positions from `octa_prowb`.
- `len_octa`: number of selected positions per sequence.

## Preprocessing

Training starts in `sitega_train()`.

1. Foreground and background FASTA files are scanned to count valid sequences.
2. Sequence lengths are loaded, filtering out sequences shorter than
   `motif_len` or longer than `max_peak_len`.
3. Foreground sequences are encoded as dinucleotide arrays for both strands.
4. Background sequences are used to estimate:
   - dinucleotide frequency means for every LPD length;
   - diagonal covariance terms for every LPD length;
   - expected dinucleotide frequencies;
   - k-mer frequencies for the `olig_bg` background model.
5. Foreground and background k-mer frequencies are converted to log-ratios:
   `log10(foreground_frequency) - log10(background_frequency)`.
6. For each possible motif window, the implementation averages the k-mer
   log-ratios across the window. This produces `octa_prow`.
7. For each sequence, positions in roughly the top third of `octa_prow` are
   stored as candidate motif positions.

Implementation notes:

- Foreground/background k-mer frequencies are counted on both strands, but
  `octa_pro1` and `octa_prow` are stored in forward-sequence coordinates.
- The current code contains a bug while computing the maximum positive position
  weight. It assigns `dw = maxw` instead of `maxw = dw`, so `maxw` remains zero
  and the intended rescaling of better positions is disabled.
- Candidate positions are sampled from the forward-coordinate high-scoring list.
  If `ori == 1`, feature scoring uses the same stored `pos` on the reverse
  complement, while `fpr` mirrors the position back to forward coordinates.
  This strand/coordinate behavior should be preserved only deliberately.

## Individual Initialization

Each individual is initialized by `town::init_rand_hoxa()`.

Per sequence:

1. A strand is selected randomly.
2. A motif position is sampled from the high-scoring k-mer position list
   `octa_prowb`, using cumulative weights `octa_prows`.

For LPD features:

1. The individual contains `size` LPDs.
2. Dinucleotide codes are sampled randomly, with a per-code cap of
   `motif_len - 1`.
3. For each selected dinucleotide code, random local positions in the motif are
   chosen.
4. Single-position LPDs may be expanded left or right up to `max_lpd`.
5. The feature list is kept sorted by dinucleotide code and start position.
6. `deg[]` and `odg[]` are updated to match the generated features.

The initialized candidate is rejected if:

- its internal consistency check fails;
- it duplicates an already initialized population member.

Accepted candidates are scored with `EvalMahFIT()` and then copied into
`pop[]`. After initialization, the population is sorted by descending fitness.

## Fitness Function

Fitness is computed by `EvalMahFIT()`.

For one candidate model and one aligned sequence window, each LPD feature
produces a value:

```text
feature_value = matching_dinucleotide_count / LPD_length
```

For all training sequences:

1. The implementation computes the foreground mean feature vector `f1`.
2. It computes the foreground feature covariance from outer products of these
   feature vectors.
3. It subtracts `f1 * f1^T` from the accumulated second moments.
4. It adds background diagonal covariance terms from `dcv`.
5. It averages foreground and background covariance by dividing by `2`.
6. It computes the difference vector:

```text
df = foreground_mean - background_mean
```

where background means come from `dav`.

If a selected LPD interval contains an invalid dinucleotide code (`-1`) in any
training sequence, the current implementation sets the candidate fitness to
zero.

The main separation score is Mahalanobis-like:

```text
mah = df^T * inverse(covariance) * df
```

The implementation explicitly inverts the covariance matrix through `BackMat()`
and then multiplies by `df`.

The position/background-bias factor is:

```text
fpr = 10 ^ mean(octa_prow[selected_position])
```

For reverse-strand placements, the position is converted back to the forward
coordinate system before reading `octa_prow`.

The final fitness used for selection is:

```text
fit = mah * fpr
```

Implementation note: `infc` controls computation of an information-content
term `inf`, but the current final fitness does not multiply by `inf`. As
implemented, `inf` is computed but not used in selection.

## Mutation Phase

Each generation starts with a mutation phase.

For each population member:

1. The current individual is copied into a temporary candidate.
2. One mutation type is selected, depending on per-individual stop flags.
3. The mutation is applied.
4. The affected LPD is re-ordered if needed.
5. The candidate is rejected if it duplicates another population member.
6. The candidate is scored with `EvalMahFIT()`.
7. The candidate replaces the parent only if its fitness is strictly higher.

Mutation types:

- `MutOlig0`: changes the dinucleotide code of one LPD while avoiding overlap
  conflicts with existing LPDs of the new dinucleotide type.
- `MutCry0`: changes the local interval of one LPD, preserving non-overlap
  constraints with neighboring LPDs of the same dinucleotide type.
- `MutRegShiftHoxaW`: changes one sequence's motif position and possibly
  strand, sampling the new position from high-scoring k-mer windows.

The mutation loop tracks successful and attempted mutations per mutation type.
If the recent success ratio for a mutation type drops below a threshold, that
mutation type is stopped for the current individual. The phase also stops when
the fitness growth per cycle becomes too small or the configured attempt budget
is reached.

This is an elitist hill-climbing mutation step: neutral and worse candidates are
discarded immediately.

## Recombination Phase

After mutation, the population is sorted again and recombination starts.

The implementation builds a weighted directed pair list from population ranks.
Better-ranked individuals receive more pairing opportunities. The pair list is
shuffled each recombination cycle.

For each selected pair:

1. Both parents are copied into temporary candidates.
2. One of six recombination operators is selected randomly.
3. The operator modifies one or both candidates.
4. Affected features are re-ordered when needed.
5. Candidates that duplicate existing population members are rejected.
6. Both candidates are scored with `EvalMahFIT()`.
7. Each candidate independently replaces its corresponding parent only if its
   fitness is strictly higher.

Recombination operators:

- `Reco2_Original`: swaps one compatible LPD between two candidates.
- `Reco2_Economic`: swaps a similar but non-identical LPD of the same
  dinucleotide type.
- `Reco2_One_dinucleotide_local`: recombines local intervals within one
  dinucleotide type.
- `Reco2_One_dinucleotide_full`: swaps all LPDs of one dinucleotide type when
  counts match.
- `Reco2_Two_dinucleotides_full`: swaps the full feature groups for two
  dinucleotide types, choosing the pair of types so each candidate keeps the
  same total number of features.
- `Reco2Peak`: swaps selected sequence positions and strands between two
  candidates.

The recombination phase stops when there are no useful local recombinations,
success ratios fall below thresholds, recombination no longer improves the best
fitness enough, or the recombination budget is exhausted.

## Generation Loop and Selection

The high-level loop is:

1. Initialize a population at the start of the run.
2. Sort by descending fitness.
3. Run mutation on each individual.
4. Sort by descending fitness.
5. Run recombination over weighted rank pairs.
6. Sort by descending fitness.
7. Measure generation-level fitness improvement.
8. Stop when improvement falls below the exit threshold after generation 2.

Selection is always elitist and local:

- candidates replace parents only when `candidate.fit > parent.fit`;
- population order is maintained by sorting after mutation/recombination;
- there is no probabilistic acceptance of worse candidates;
- there is no explicit crossover generation replacement step.

The best individual is therefore monotonically non-decreasing in fitness within
a run, except for possible implementation bugs or numerical failures.

## Final Output

After the generation loop exits, the implementation writes model files for
every population member:

- `*_matN`: LPD feature matrix and feature weights;
- `*_locN`: best-scoring motif location per sequence for that model.

Final feature weights are computed by `EvalMahFITTrain()`, which recomputes the
Mahalanobis model and scans every possible motif window on both strands for each
sequence.

Implementation notes:

- The Python wrapper usually consumes only the first few motifs, but the C++
  code currently writes files for all `pop_size` individuals.
- `EvalMahFITTrain()` appears to omit normalization of the foreground feature
  sums before covariance and `df` computation. The final matrix/location output
  may therefore be mathematically inconsistent with the training-time
  `EvalMahFIT()` score.

## Mathematical Summary

The model searches for a set of local dinucleotide features that separates
foreground motif-aligned windows from background dinucleotide statistics.

The key ideas are:

- represent a motif by multiple local dinucleotide-presence features rather
  than by a standard position weight matrix;
- assign each training sequence a current motif window and strand;
- evaluate a candidate by how far its foreground feature mean is from the
  background mean under a covariance-normalized distance;
- bias motif window placement toward regions enriched for foreground k-mers;
- optimize the model with greedy genetic operations over LPD definitions and
  per-sequence placements.

The core score is:

```text
fit = (df^T * covariance^-1 * df) * 10^mean(kmer_log_ratio)
```

where:

- `df` is the foreground-background difference in LPD feature means;
- `covariance` is the averaged foreground/background covariance estimate;
- `mean(kmer_log_ratio)` rewards motif placements in foreground-enriched
  sequence windows.

## Current Implementation Caveats

The implementation is old C-style C++ and should be treated carefully when
modifying:

- it relies on fixed-size global arrays and manual memory management;
- many operations copy full `pos[]` and `ori[]` arrays;
- fitness is recomputed from scratch for most candidate changes;
- `BackMat()` computes a full matrix inverse although only a linear solve is
  mathematically required;
- `infc` currently does not affect the final `fit`;
- the weighted position initialization contains a `maxw` assignment bug;
- some recombination duplicate checks use stale pair indices after pair-list
  shuffling;
- `Reco2_Two_dinucleotides_full()` appears to update `odg[]` using positions at
  the beginning of `tot[]`, not the selected dinucleotide groups.
