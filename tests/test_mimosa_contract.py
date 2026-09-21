from __future__ import annotations

import json

import numpy as np
import pytest
from mimosa import (
    BestPerSequence,
    EncodedSequences,
    TopFractionHits,
    read_fasta,
    read_model,
    reconstruct_pfm,
    scan,
    write_model,
)


def test_fasta_names_and_order_are_preserved(tmp_path):
    path = tmp_path / "input.fa"
    path.write_text(">first description\nACGT\n>second\nNNAA\n")

    sequences, names = read_fasta(path)

    assert names == ("first", "second")
    subset = EncodedSequences.from_rows(sequences[index] for index in [1, 0])
    assert subset[0].tolist() == [4, 4, 0, 0]
    assert subset[1].tolist() == [0, 1, 2, 3]


def test_scan_and_pfm_reconstruction_use_mimosa(pwm_model, sequences):
    scores = scan(pwm_model, sequences, strands="best")
    pfm = reconstruct_pfm(
        pwm_model,
        sequences,
        TopFractionHits(0.10, BestPerSequence()),
        strands="best",
    )

    assert len(scores) == len(sequences)
    assert scores.data.dtype == np.float32
    assert pfm.shape == (4, pwm_model.motif_length)
    np.testing.assert_allclose(pfm.sum(axis=0), 1.0, atol=1e-6)


def test_model_bundle_round_trip_and_no_overwrite(tmp_path, pwm_model):
    path = tmp_path / "model"
    write_model(path, pwm_model)
    loaded = read_model(path)

    assert loaded.name == pwm_model.name
    np.testing.assert_array_equal(loaded.weights, pwm_model.weights)
    with pytest.raises(Exception):
        write_model(path, pwm_model)
    manifest = json.loads(json.dumps({"path": str(path)}))
    assert manifest["path"] == str(path)
