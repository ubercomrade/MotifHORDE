from __future__ import annotations

import os

import numpy as np
import pytest

from hordemotifs.batches import make_sequence_batch
from hordemotifs.functions import pfm_to_pwm
from hordemotifs.models import GenericModel


@pytest.fixture
def test_pfm():
    return np.array(
        [
            [0.8, 0.1, 0.1, 0.1, 0.2, 0.25],
            [0.05, 0.8, 0.05, 0.1, 0.2, 0.25],
            [0.05, 0.05, 0.8, 0.1, 0.3, 0.25],
            [0.1, 0.05, 0.05, 0.7, 0.3, 0.25],
        ],
        dtype=np.float32,
    )


@pytest.fixture
def pwm_model(test_pfm):
    pwm = pfm_to_pwm(test_pfm)
    representation = np.concatenate((pwm, np.min(pwm, axis=0, keepdims=True)), axis=0)
    return GenericModel("pwm", "M1", representation.astype(np.float32), test_pfm.shape[1], {"kmer": 1, "_source_pfm": test_pfm})


@pytest.fixture
def sequence_batch():
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
    seqs = [
        "ACGTACGTACGT",
        "TTTTTTTTTTTT",
        "AAAAAACCCCCC",
    ]
    return make_sequence_batch(np.array([mapping[base] for base in seq], dtype=np.int8) for seq in seqs)


@pytest.fixture
def sample_meme(tmp_path, test_pfm):
    path = tmp_path / "test.meme"
    with open(path, "w") as handle:
        handle.write("MEME version 4\nALPHABET= ACGT\nstrands: + -\n\n")
        handle.write("MOTIF TestMotif1\n")
        handle.write("letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0\n")
        np.savetxt(handle, test_pfm.T, fmt="%.6f")
    return os.fspath(path)
