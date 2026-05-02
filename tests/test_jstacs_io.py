from __future__ import annotations

from motifhorde.io import write_jstacs_fasta


def test_write_jstacs_fasta_empty(tmp_path):
    source = tmp_path / "empty.fa"
    target = tmp_path / "annot.fa"
    source.write_text("")

    write_jstacs_fasta(source, target)

    assert target.read_text() == ""


def test_write_jstacs_fasta_single_sequence(tmp_path):
    source = tmp_path / "input.fa"
    target = tmp_path / "annot.fa"
    source.write_text(">x\nACGTAC\n")

    write_jstacs_fasta(source, target)

    assert target.read_text() == "> position: 3; value: 1.0\nACGTAC\n"


def test_write_jstacs_fasta_multiple_sequences(tmp_path):
    source = tmp_path / "input.fa"
    target = tmp_path / "annot.fa"
    source.write_text(">x\nACGT\n>y\nTTTTAA\n")

    write_jstacs_fasta(source, target)

    assert target.read_text() == (
        "> position: 2; value: 1.0\n"
        "ACGT\n"
        "> position: 3; value: 2.0\n"
        "TTTTAA\n"
    )


def test_write_jstacs_fasta_custom_tags(tmp_path):
    source = tmp_path / "input.fa"
    target = tmp_path / "annot.fa"
    source.write_text(">x\nACGT\n")

    write_jstacs_fasta(source, target, position_tag="peak", value_tag="signal")

    assert target.read_text() == "> peak: 2; signal: 1.0\nACGT\n"
