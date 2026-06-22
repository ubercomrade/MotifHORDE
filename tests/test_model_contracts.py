from __future__ import annotations

import os

import joblib
import numpy as np
import pytest
from mimosa.batches import make_sequence_batch

from motifhorde.comparison import TomtomComparator, UniversalMotifComparator
from motifhorde.models import (
    GenericModel,
    calculate_threshold_table,
    get_pfm,
    get_sites,
    read_model,
    scan_model,
    scan_model_strands,
    write_model,
)
from motifhorde.pipeline import DeNovoPipeline


EXPECTED_PLUS_SCORES = np.array(
    [
        [4.297346, -5.959855, -4.862407, -4.170259, 4.297346, -5.959855, -4.862407],
        [-2.917354, -2.917354, -2.917354, -2.917354, -2.917354, -2.917354, -2.917354],
        [-1.804842, -1.804842, -1.804842, -1.804842, -2.496990, -0.418423, -3.189139],
    ],
    dtype=np.float32,
)
EXPECTED_MINUS_SCORES = np.array(
    [
        [-4.862407, -5.959855, 4.297346, -4.170259, -4.862407, -5.959855, 4.297346],
        [-1.804842, -1.804842, -1.804842, -1.804842, -1.804842, -1.804842, -1.804842],
        [-2.917354, -3.609503, -3.609503, -0.838787, -2.783840, -2.783840, -2.783840],
    ],
    dtype=np.float32,
)
EXPECTED_BEST_SCORES = np.array(
    [
        [4.297346, -5.959855, 4.297346, -4.170259, 4.297346, -5.959855, 4.297346],
        [-1.804842, -1.804842, -1.804842, -1.804842, -1.804842, -1.804842, -1.804842],
        [-1.804842, -1.804842, -1.804842, -0.838787, -2.496990, -0.418423, -2.783840],
    ],
    dtype=np.float32,
)
EXPECTED_THRESHOLD_TABLE = np.array(
    [
        [4.297346, 1.021189],
        [-0.418423, 0.924279],
        [-0.838787, 0.845098],
        [-1.804842, 0.392800],
        [-2.496990, 0.367977],
        [-2.783840, 0.301030],
        [-2.917354, 0.160851],
        [-3.189139, 0.146128],
        [-3.609503, 0.118099],
        [-4.170259, 0.091770],
        [-4.862407, 0.043466],
        [-5.959855, 0.0],
    ],
    dtype=np.float32,
)
EXPECTED_BEST_PFM = np.array(
    [
        [0.8125, 0.3125, 0.3125, 0.3125, 0.5625, 0.3125],
        [0.0625, 0.5625, 0.3125, 0.3125, 0.3125, 0.5625],
        [0.0625, 0.0625, 0.3125, 0.0625, 0.0625, 0.0625],
        [0.0625, 0.0625, 0.0625, 0.3125, 0.0625, 0.0625],
    ],
    dtype=np.float32,
)
EXPECTED_TOP_FRACTION_PFM = np.array(
    [
        [0.625, 0.125, 0.125, 0.125, 0.625, 0.125],
        [0.125, 0.625, 0.125, 0.125, 0.125, 0.625],
        [0.125, 0.125, 0.625, 0.125, 0.125, 0.125],
        [0.125, 0.125, 0.125, 0.625, 0.125, 0.125],
    ],
    dtype=np.float32,
)


def _copy_model(
    model: GenericModel, name: str, *, type_key: str | None = None
) -> GenericModel:
    representation = (
        None
        if model.representation is None
        else np.array(model.representation, copy=True)
    )
    config = {
        key: np.array(value, copy=True) if isinstance(value, np.ndarray) else value
        for key, value in model.config.items()
    }
    return GenericModel(
        type_key or model.type_key, name, representation, model.length, config
    )


def _write_bamm(path) -> None:
    path.write_text(
        "\n\n".join(
            [
                "0.70 0.10 0.10 0.10",
                "0.10 0.70 0.10 0.10",
                "0.10 0.10 0.70 0.10",
                "0.10 0.10 0.10 0.70",
            ]
        )
    )


def _write_sitega(path) -> None:
    path.write_text(
        "\n".join(
            [
                "SiteGAContract",
                "1\tLPD count",
                "4\tModel length",
                "0.0\tMinimum",
                "1.0\tRazmah",
                "0\t2\t1.5\t0\tac",
            ]
        )
    )


def _write_scores(path) -> None:
    path.write_text(">seq0\n1.0 2.0 3.0\n>seq1\n0.5 0.25\n")


def test_mimosa_dependency_contract():
    import mimosa
    from mimosa.models import write_model as mimosa_write_model
    from mimosa.scanning import calculate_threshold_table as mimosa_threshold_table
    from mimosa.scanning import scan_model_strands as mimosa_scan_model_strands

    assert GenericModel is mimosa.GenericModel
    assert callable(mimosa.scan_model)
    assert callable(mimosa.get_sites)
    assert callable(mimosa.get_pfm)
    assert callable(mimosa_write_model)
    assert callable(mimosa_threshold_table)
    assert callable(mimosa_scan_model_strands)


def test_scan_model_contract_values(pwm_model, sequence_batch):
    expected_by_strand = {
        "+": EXPECTED_PLUS_SCORES,
        "-": EXPECTED_MINUS_SCORES,
        "best": EXPECTED_BEST_SCORES,
    }

    for strand, expected_values in expected_by_strand.items():
        scores = scan_model(pwm_model, sequence_batch, strand=strand)

        assert scores["values"].dtype == np.float32
        assert scores["mask"].dtype == bool
        assert scores["lengths"].tolist() == [7, 7, 7]
        assert scores["mask"].all()
        np.testing.assert_allclose(
            scores["values"], expected_values, rtol=1e-6, atol=1e-6
        )

    strand_bundle = scan_model_strands(pwm_model, sequence_batch)
    assert strand_bundle["values"].shape == (2, 3, 7)
    assert strand_bundle["lengths"].tolist() == [7, 7, 7]
    np.testing.assert_allclose(
        strand_bundle["values"][0], EXPECTED_PLUS_SCORES, rtol=1e-6, atol=1e-6
    )
    np.testing.assert_allclose(
        strand_bundle["values"][1], EXPECTED_MINUS_SCORES, rtol=1e-6, atol=1e-6
    )


def test_site_table_contract_schema_and_values(pwm_model, sequence_batch):
    sites = get_sites(pwm_model, sequence_batch)

    assert list(sites.columns) == [
        "seq_index",
        "start",
        "end",
        "strand",
        "score",
        "log_tail",
        "site",
    ]
    assert sites["seq_index"].tolist() == [0, 1, 2]
    assert sites["start"].tolist() == [0, 0, 5]
    assert sites["end"].tolist() == [6, 6, 11]
    assert sites["strand"].tolist() == ["+", "-", "+"]
    assert sites["site"].tolist() == ["ACGTAC", "AAAAAA", "ACCCCC"]
    np.testing.assert_allclose(
        sites["score"], [4.297346, -1.804842, -0.418423], rtol=1e-6, atol=1e-6
    )
    np.testing.assert_allclose(
        sites["log_tail"], [1.021189, 0.392800, 0.924279], rtol=1e-6, atol=1e-6
    )


def test_threshold_and_pfm_contracts(pwm_model, sequence_batch):
    threshold_table = calculate_threshold_table(pwm_model, sequence_batch)
    np.testing.assert_allclose(
        threshold_table, EXPECTED_THRESHOLD_TABLE, rtol=1e-6, atol=1e-6
    )

    with pytest.raises(ValueError, match="fpr_threshold is required"):
        get_sites(pwm_model, sequence_batch, mode="threshold")

    threshold_sites = get_sites(
        pwm_model,
        sequence_batch,
        mode="threshold",
        fpr_threshold=0.5,
        threshold_table=threshold_table,
    )
    assert threshold_sites["site"].tolist() == [
        "ACGTAC",
        "ACGTAC",
        "ACGTAC",
        "ACGTAC",
        "AAAAAA",
        "AAAAAA",
        "AAAAAA",
        "AAAAAA",
        "AAAAAA",
        "AAAAAA",
        "AAAAAA",
        "ACCCCC",
        "GGGTTT",
        "AAAAAA",
        "AAAAAC",
        "AAAACC",
        "AAACCC",
        "AACCCC",
        "GGGGTT",
        "GGGGGT",
        "GGGGGG",
    ]

    pfm = get_pfm(pwm_model, sequence_batch)
    np.testing.assert_allclose(pfm, EXPECTED_BEST_PFM, rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(
        pfm.sum(axis=0), np.ones(pwm_model.length), rtol=1e-6, atol=1e-6
    )

    top_fraction_pfm = get_pfm(pwm_model, sequence_batch, top_fraction=0.5)
    np.testing.assert_allclose(
        top_fraction_pfm, EXPECTED_TOP_FRACTION_PFM, rtol=1e-6, atol=1e-6
    )


def test_empty_and_short_sequence_contracts(pwm_model):
    empty_batch = make_sequence_batch([])
    empty_scores = scan_model(pwm_model, empty_batch, strand="best")
    assert empty_scores["values"].shape == (0, 0)
    assert empty_scores["mask"].shape == (0, 0)
    assert empty_scores["lengths"].shape == (0,)

    short_batch = make_sequence_batch([np.array([0, 1, 2], dtype=np.int8)])
    short_scores = scan_model(pwm_model, short_batch, strand="best")
    assert short_scores["values"].shape == (1, 0)
    assert short_scores["mask"].shape == (1, 0)
    assert short_scores["lengths"].tolist() == [0]

    short_sites = get_sites(pwm_model, short_batch)
    assert list(short_sites.columns) == [
        "seq_index",
        "start",
        "end",
        "strand",
        "score",
        "log_tail",
        "site",
    ]
    assert short_sites.empty
    with pytest.raises(ValueError, match="No sites found"):
        get_pfm(pwm_model, short_batch)


def test_read_model_contract_for_supported_families(
    tmp_path, pwm_model, sequence_batch, sample_meme
):
    bamm_path = tmp_path / "contract.ihbcp"
    sitega_path = tmp_path / "contract.mat"
    dimont_path = tmp_path / "contract_dimont.pkl"
    slim_path = tmp_path / "contract_slim.pkl"
    scores_path = tmp_path / "contract_scores.txt"

    _write_bamm(bamm_path)
    _write_sitega(sitega_path)
    _write_scores(scores_path)
    joblib.dump(
        _copy_model(pwm_model, "DimontContract", type_key="dimont"), dimont_path
    )
    joblib.dump(_copy_model(pwm_model, "SlimContract", type_key="slim"), slim_path)

    models = {
        "pwm": read_model(sample_meme, "pwm"),
        "bamm": read_model(os.fspath(bamm_path), "bamm", order=0),
        "sitega": read_model(os.fspath(sitega_path), "sitega"),
        "dimont": read_model(os.fspath(dimont_path), "dimont"),
        "slim": read_model(os.fspath(slim_path), "slim"),
        "scores": read_model(os.fspath(scores_path), "scores"),
    }

    assert models["pwm"].name == "TestMotif1"
    assert models["bamm"].config["order"] == 0
    assert models["sitega"].config["kmer"] == 2
    assert models["dimont"].name == "DimontContract"
    assert models["slim"].name == "SlimContract"
    assert models["scores"].name == "contract_scores"

    for type_key, model in models.items():
        assert model.type_key == type_key
        assert isinstance(model.config, dict)

        scores = scan_model(
            model, None if type_key == "scores" else sequence_batch, strand="best"
        )
        assert {"values", "lengths", "padding_value"}.issubset(scores)
        if type_key != "scores":
            assert scores["values"].shape[0] == len(sequence_batch["lengths"])

        strand_bundle = scan_model_strands(
            model, None if type_key == "scores" else sequence_batch
        )
        assert strand_bundle["values"].shape[0] == 2
        assert strand_bundle["values"].shape[1] == len(strand_bundle["lengths"])


def test_read_model_rejects_legacy_pwm_txt_extension(tmp_path, sample_meme):
    txt_path = tmp_path / "legacy_pwm.txt"
    with open(sample_meme) as source:
        txt_path.write_text(source.read())

    with pytest.raises(ValueError):
        read_model(os.fspath(txt_path), "pwm")


def test_pickle_and_pwm_write_contract(tmp_path, pwm_model, test_pfm):
    pickle_path = tmp_path / "model.pkl"
    joblib.dump(pwm_model, pickle_path)

    loaded_pickle = read_model(os.fspath(pickle_path), "pwm")
    assert isinstance(loaded_pickle, GenericModel)
    assert loaded_pickle.type_key == "pwm"
    assert loaded_pickle.name == pwm_model.name
    np.testing.assert_allclose(
        loaded_pickle.config["_source_pfm"], test_pfm, rtol=1e-6, atol=1e-6
    )

    pfm_path = tmp_path / "written.pfm"
    write_model(pwm_model, os.fspath(pfm_path))
    loaded_pfm = read_model(os.fspath(pfm_path), "pwm")
    assert loaded_pfm.name == "written"
    assert loaded_pfm.length == pwm_model.length
    np.testing.assert_allclose(
        loaded_pfm.config["_source_pfm"], test_pfm, rtol=1e-6, atol=1e-6
    )


def test_read_model_rejects_unsupported_pickle_payload(tmp_path):
    pickle_path = tmp_path / "unsupported.pkl"
    joblib.dump({"type_key": "pwm"}, pickle_path)

    with pytest.raises(TypeError, match="Unsupported PWM pickle payload"):
        read_model(os.fspath(pickle_path), "pwm")


def test_pipeline_saved_models_are_mimosa_readable(
    tmp_path, pwm_model, sequence_batch, test_pfm
):
    pipeline = DeNovoPipeline(
        discovery_tool=None,
        evaluator=None,
        comparator=TomtomComparator(n_jobs=1),
    )
    final_info = [(pwm_model.name, {"length": pwm_model.length})]
    final_stats = {
        f"{pwm_model.name}_length-{pwm_model.length}": {
            "auPRC": 0.7,
            "auROC": 0.8,
            "pauPRC": 0.6,
            "pauROC": 0.9,
        }
    }

    pipeline._save_results(
        [pwm_model],
        final_info,
        final_stats,
        os.fspath(tmp_path),
        "pauROC",
        sequence_batch,
    )

    model_path = tmp_path / "models" / f"001_{pwm_model.name}.pkl"
    loaded_model = read_model(os.fspath(model_path), "pwm")
    assert isinstance(loaded_model, GenericModel)
    assert loaded_model.type_key == "pwm"
    assert loaded_model.name == pwm_model.name
    np.testing.assert_allclose(
        loaded_model.config["_source_pfm"], test_pfm, rtol=1e-6, atol=1e-6
    )

    meme_path = tmp_path / "models" / "all_motifs_in_pfm_form.meme"
    loaded_meme = read_model(os.fspath(meme_path), "pwm")
    assert loaded_meme.type_key == "pwm"
    assert loaded_meme.name == pwm_model.name
    np.testing.assert_allclose(
        loaded_meme.config["_source_pfm"], test_pfm, rtol=1e-6, atol=1e-6
    )


def test_motif_comparison_contract_same_type_and_ordering(pwm_model):
    copy_model = _copy_model(pwm_model, "M2")
    shifted_model = _copy_model(pwm_model, "M3")
    shifted_model.representation = np.roll(shifted_model.representation, 1, axis=1)
    shifted_model.config["_source_pfm"] = np.roll(
        shifted_model.config["_source_pfm"], 1, axis=1
    )

    frame = TomtomComparator(n_jobs=1).compare(
        [pwm_model],
        [copy_model, shifted_model],
    )

    assert frame["query"].tolist() == ["M1", "M1"]
    assert frame["target"].tolist() == ["M2", "M3"]
    assert frame["orientation"].tolist() == ["+-", "++"]
    assert frame["offset"].tolist() == [-2, -1]
    np.testing.assert_allclose(frame["score"], [0.937834, 1.0], rtol=1e-6, atol=1e-6)


def test_motif_comparison_contract_heterogeneous_pfm_reconstruction(
    pwm_model, sequence_batch
):
    heterogeneous_model = _copy_model(pwm_model, "H1", type_key="dimont")

    frame = TomtomComparator(n_jobs=1).compare(
        [pwm_model],
        [heterogeneous_model],
        sequences=sequence_batch,
    )

    assert frame.loc[0, "query"] == "M1"
    assert frame.loc[0, "target"] == "H1"
    assert frame.loc[0, "orientation"] == "++"
    assert frame.loc[0, "offset"] == 0
    assert frame.loc[0, "metric"] == "pcc"
    assert frame.loc[0, "score"] == pytest.approx(1.0)


def test_profile_comparison_contract_fixed_seed(pwm_model, sequence_batch):
    target_model = _copy_model(pwm_model, "M2")
    comparator = UniversalMotifComparator(
        metric="co",
        n_jobs=1,
        seed=123,
    )

    frame = comparator.compare([pwm_model], [target_model], sequences=sequence_batch)

    assert frame.loc[0, "query"] == "M1"
    assert frame.loc[0, "target"] == "M2"
    assert frame.loc[0, "orientation"] == "++"
    assert frame.loc[0, "offset"] == 0
    assert frame.loc[0, "metric"] == "co"
    assert frame.loc[0, "n_sites"] == 0
    assert frame.loc[0, "score"] == pytest.approx(0.0)
