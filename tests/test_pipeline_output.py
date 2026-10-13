from __future__ import annotations

import json
import os
from dataclasses import replace

import numpy as np
import pytest
from mimosa import PWM, EncodedSequences, build_null, read_model
from mimosa.io import write_null_bundle

from motifhorde.comparison import ComparisonSettings, comparison_sequences
from motifhorde.discovery import MotifDiscoveryTool
from motifhorde.evaluation import EvaluatedMotif, PerformanceEvaluator
from motifhorde.pipeline import DeNovoPipeline, _select_nonredundant


class FakeDiscovery(MotifDiscoveryTool):
    def __init__(self, model):
        super().__init__("fake")
        self.model = model

    def discover(
        self, foreground, background, output_dir, number_of_motifs, *args, **kwargs
    ):
        return [
            replace(self.model, name=f"Fake-{index}")
            for index in range(1, number_of_motifs + 1)
        ]


def _sequences_with_rare_sites(site: str) -> tuple[PWM, EncodedSequences]:
    # Mimosa keeps only positions with -log10(ERR) >= 4, and ERR is empirical
    # on the compared sequences: a site must be rarer than 1 in 10,000
    # positions, otherwise no sites remain and every score collapses to zero.
    weights = np.full((5, len(site)), -2.0, dtype=np.float32)
    weights[4] = -3.0
    for position, base in enumerate(site):
        weights["ACGT".index(base), position] = 3.0
    rng = np.random.default_rng(5)
    rows = []
    for index in range(400):
        row = list(rng.choice(list("ACGT"), 200))
        if index < 4:
            start = int(rng.integers(0, 200 - len(site)))
            row[start : start + len(site)] = site
        rows.append("".join(row))
    return PWM("Strong", weights, (0.25,) * 4), EncodedSequences.from_strings(rows)


def test_select_nonredundant_keeps_best_of_similar_motifs():
    pwm_model, long_sequences = _sequences_with_rare_sites("ACGTTGCA")
    worse = EvaluatedMotif(
        replace(pwm_model, name="Similar"),
        "worse",
        "task_0000",
        "odd",
        {},
        {"pauROC": 0.4},
    )
    best = EvaluatedMotif(pwm_model, "best", "task_0001", "odd", {}, {"pauROC": 0.9})
    settings = ComparisonSettings(criterion="score", threshold=0.9)

    kept = _select_nonredundant(
        [
            (worse.model, worse.metrics["pauROC"], worse),
            (best.model, best.metrics["pauROC"], best),
        ],
        long_sequences,
        settings,
    )

    assert [record.id for record in kept] == ["best"]


def test_pipeline_fills_default_comparison_threshold():
    pipeline = DeNovoPipeline(FakeDiscovery(None), PerformanceEvaluator("peaks"))

    assert pipeline.comparison.threshold == pytest.approx(0.9)


def test_pipeline_writes_versioned_bundles_and_refuses_overwrite(tmp_path, pwm_model):
    foreground = tmp_path / "foreground.fa"
    background = tmp_path / "background.fa"
    promoters = tmp_path / "promoters.fa"
    for path in (foreground, background, promoters):
        path.write_text(
            ">one\nACGTACGTACGT\n>two\nTACGTACGTACG\n>three\nACGTACGTACGT\n>four\nTACGTACGTACG\n"
        )
    output = tmp_path / "output"
    pipeline = DeNovoPipeline(
        FakeDiscovery(pwm_model),
        PerformanceEvaluator("peaks"),
        comparison_threshold=-1,
        number_of_motifs=1,
    )

    pipeline.run(
        os.fspath(foreground),
        os.fspath(background),
        os.fspath(promoters),
        os.fspath(output),
        {"length": [4]},
        "pauROC",
    )

    root = output / "fake"
    bootstrap_stats = json.loads((root / "bootstrap" / "statistics.json").read_text())
    final_stats = json.loads((root / "motifs" / "statistics.json").read_text())
    assert bootstrap_stats["schema_version"] == 1
    assert bootstrap_stats["motifs"][0]["split"] == "odd"
    assert final_stats["schema_version"] == 1
    assert list((root / "bootstrap" / "models").iterdir())
    assert not list(root.rglob("*.pkl"))
    for entry in bootstrap_stats["motifs"]:
        assert read_model(root / "bootstrap" / entry["model"]).name == entry["name"]

    with pytest.raises(FileExistsError):
        pipeline.run(
            os.fspath(foreground),
            os.fspath(background),
            os.fspath(promoters),
            os.fspath(output),
            {"length": [4]},
            "pauROC",
        )


def test_pipeline_seed_and_order_are_independent_of_jobs(tmp_path, pwm_model):
    foreground = tmp_path / "foreground.fa"
    background = tmp_path / "background.fa"
    promoters = tmp_path / "promoters.fa"
    for path in (foreground, background, promoters):
        path.write_text(
            ">one\nACGTACGTACGT\n>two\nTACGTACGTACG\n>three\nACGTACGTACGT\n>four\nTACGTACGTACG\n"
        )

    def run(output, jobs):
        DeNovoPipeline(
            FakeDiscovery(pwm_model),
            PerformanceEvaluator("peaks"),
            comparison_threshold=-1,
            number_of_motifs=1,
            jobs=jobs,
            seed=41,
        ).run(
            os.fspath(foreground),
            os.fspath(background),
            os.fspath(promoters),
            os.fspath(output),
            {"length": [4]},
            "pauROC",
        )

    run(tmp_path / "jobs1", 1)
    run(tmp_path / "jobs2", 2)
    for relative in (
        "fake/bootstrap/statistics.json",
        "fake/motifs/statistics.json",
    ):
        first = (tmp_path / "jobs1" / relative).read_bytes()
        second = (tmp_path / "jobs2" / relative).read_bytes()
        assert first == second


def test_pipeline_accepts_mimosa_pvalue_bundle(tmp_path, pwm_model):
    foreground = tmp_path / "foreground.fa"
    background = tmp_path / "background.fa"
    promoters = tmp_path / "promoters.fa"
    fasta = ">one\nACGTACGTACGT\n>two\nTACGTACGTACG\n>three\nACGTACGTACGT\n>four\nTACGTACGTACG\n"
    for path in (foreground, background, promoters):
        path.write_text(fasta)

    other = replace(pwm_model, name="Other")
    null_path = tmp_path / "null"
    write_null_bundle(
        null_path,
        build_null(
            [pwm_model, other],
            sequences=comparison_sequences(),
            n_samples=8,
            seed=13,
            search_range=10,
        ),
    )

    DeNovoPipeline(
        FakeDiscovery(pwm_model),
        PerformanceEvaluator("peaks"),
        comparison_criterion="p-value",
        comparison_threshold=1.0,
        null_distribution=os.fspath(null_path),
        number_of_motifs=1,
    ).run(
        os.fspath(foreground),
        os.fspath(background),
        os.fspath(promoters),
        os.fspath(tmp_path / "output"),
        {"length": [4]},
        "pauROC",
    )

    assert (tmp_path / "output/fake/motifs/statistics.json").is_file()
