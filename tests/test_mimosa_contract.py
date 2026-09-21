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


def test_all_discovery_readers_return_mimosa_models(tmp_path):
    files = {
        "meme": (
            "motif.meme",
            "MEME version 4\n\nALPHABET= ACGT\n\n"
            "MOTIF Fixture\n"
            "letter-probability matrix: alength= 4 w= 4 nsites= 1 E= 0\n"
            + "0.25 0.25 0.25 0.25\n"
            * 4,
        ),
        "pfm": ("motif.pfm", "0.25 0.25 0.25 0.25\n" * 4),
        "bamm": (
            "motif.ihbcp",
            "\n\n".join(["0.25 0.25 0.25 0.25"] * 4),
        ),
        "sitega": (
            "motif.mat",
            "SiteGATest\n1\tLPD count\n4\tModel length\n"
            "0.0\tMinimum\n1.0\tRazmah\n0\t2\t1.5\t0\tac\n",
        ),
        "dimont": (
            "motif.xml",
            "<root><MarkovModelDiffSM><bayesianNetworkSF><trees>"
            "<pos><parameterTree><root><treeElement><contextPos>0</contextPos>"
            "<pars><pos><parameter><symbol>0</symbol><value>0.25</value></parameter>"
            "</pos><pos><parameter><symbol>1</symbol><value>0.25</value></parameter>"
            "</pos><pos><parameter><symbol>2</symbol><value>0.25</value></parameter>"
            "</pos><pos><parameter><symbol>3</symbol><value>0.25</value></parameter>"
            "</pos></pars></treeElement></root></parameterTree></pos>"
            "</trees></bayesianNetworkSF></MarkovModelDiffSM></root>",
        ),
        "slim": (
            "slim.xml",
            "<root><SLIM><length>1</length>"
            "<componentMixtureParameters><pos><pos>1</pos></pos>"
            "</componentMixtureParameters>"
            "<ancestorMixtureParameters><pos><pos><pos>1</pos></pos></pos>"
            "</ancestorMixtureParameters>"
            "<dependencyParameters><pos><pos><pos><pos>0.25</pos><pos>0.25</pos>"
            "<pos>0.25</pos><pos>0.25</pos></pos></pos></pos>"
            "</dependencyParameters></SLIM></root>",
        ),
    }

    expected_lengths = {
        "meme": 4,
        "pfm": 4,
        "bamm": 4,
        "sitega": 4,
        "dimont": 1,
        "slim": 1,
    }
    for model_format, (filename, content) in files.items():
        path = tmp_path / filename
        path.write_text(content)
        kwargs = {"format": model_format}
        if model_format == "bamm":
            kwargs["order"] = 0
        model = read_model(path, **kwargs)
        assert type(model).__name__ in {"PWM", "BaMM", "SiteGA", "Dimont", "Slim"}
        assert model.motif_length == expected_lengths[model_format]
