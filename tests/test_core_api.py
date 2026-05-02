from __future__ import annotations

import numpy as np

from motifhorde.batches import make_sequence_batch
from motifhorde.comparison import TomtomComparator, UniversalMotifComparator
from motifhorde.evaluation import PerformanceEvaluator
from motifhorde.functions import (
    batch_all_scores,
    batch_all_scores_strands,
    build_score_log_tail_table,
    lookup_score_for_tail_probability,
)
from motifhorde.io import read_fasta, write_fasta
from motifhorde.models import (
    GenericModel,
    calculate_threshold_table,
    get_frequencies,
    get_pfm,
    get_sites,
    read_model,
    scan_model,
    scan_model_strands,
)


def test_make_sequence_batch_empty():
    batch = make_sequence_batch([])
    assert batch["values"].shape == (0, 0)
    assert batch["lengths"].shape == (0,)


def test_read_write_fasta_roundtrip(tmp_path, sequence_batch):
    path = tmp_path / "seqs.fa"
    write_fasta(sequence_batch, path)
    observed = read_fasta(path)
    np.testing.assert_array_equal(observed["values"], sequence_batch["values"])
    np.testing.assert_array_equal(observed["lengths"], sequence_batch["lengths"])


def test_read_model_and_scan(sample_meme, sequence_batch):
    model = read_model(sample_meme, "pwm")
    assert isinstance(model, GenericModel)
    scores = scan_model(model, sequence_batch, strand="best")
    assert scores["values"].shape[0] == len(sequence_batch["lengths"])
    assert np.any(scores["mask"])


def test_read_meme_accepts_compact_width_header(tmp_path, test_pfm):
    path = tmp_path / "compact.meme"
    with open(path, "w") as handle:
        handle.write("MEME version 4\n\n")
        handle.write("MOTIF Compact\n")
        handle.write("letter-probability matrix: alength= 4 w=6 nsites= 20 E= 0\n")
        np.savetxt(handle, test_pfm.T, fmt="%.6f")

    model = read_model(path, "pwm")
    assert model.length == test_pfm.shape[1]


def test_batch_scoring_kernels(pwm_model, sequence_batch):
    scores = batch_all_scores(sequence_batch, pwm_model.representation, kmer=1)
    plus, minus = batch_all_scores_strands(sequence_batch, pwm_model.representation, kmer=1)
    assert scores["values"].shape == plus["values"].shape
    assert plus["values"].shape == minus["values"].shape


def test_model_functional_helpers(pwm_model, sequence_batch):
    strand_bundle = scan_model_strands(pwm_model, sequence_batch)
    assert strand_bundle["values"].shape[0] == 2

    frequencies = get_frequencies(pwm_model, sequence_batch)
    assert frequencies["values"].shape[0] == len(sequence_batch["lengths"])

    sites = get_sites(pwm_model, sequence_batch)
    assert {"seq_index", "start", "strand", "score", "site"}.issubset(sites.columns)

    pfm = get_pfm(pwm_model, sequence_batch)
    assert pfm.shape == (4, pwm_model.length)

    table = calculate_threshold_table(pwm_model, sequence_batch)
    assert table.shape[1] == 2


def test_evaluator_stores_statistics(pwm_model, sequence_batch):
    evaluator = PerformanceEvaluator(background_type="peaks")
    stats = evaluator.evaluate(pwm_model, sequence_batch, sequence_batch, err_threshold=0.1)
    assert "auROC" in stats
    assert pwm_model.config["statistics"]["auROC"] == stats["auROC"]


def test_comparator_wrappers(pwm_model, sequence_batch):
    second = GenericModel(
        "pwm",
        "M2",
        pwm_model.representation.copy(),
        pwm_model.length,
        dict(pwm_model.config),
    )
    tomtom = TomtomComparator(n_permutations=0)
    motif_frame = tomtom.compare([pwm_model], [second], sequences=sequence_batch)
    assert motif_frame.loc[0, "query"] == "M1"

    continuous = UniversalMotifComparator(metric="co", n_permutations=0)
    profile_frame = continuous.compare([pwm_model], [second], sequences=sequence_batch)
    assert profile_frame.loc[0, "target"] == "M2"


def test_score_tail_table_handles_empty_input():
    table = build_score_log_tail_table(np.array([], dtype=np.float32))
    assert table.shape == (1, 2)


def test_lookup_score_for_strict_tail_probability_uses_strictest_cutoff():
    table = build_score_log_tail_table(np.array([10.0, 9.0, 8.0, 7.0], dtype=np.float32))
    score = lookup_score_for_tail_probability(table, 0.01)
    assert score == 10.0
