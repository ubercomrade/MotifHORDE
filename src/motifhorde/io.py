"""Local I/O missing from Mimosa: discovery FASTA variants and MEME export."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
from mimosa import EncodedSequences, read_fasta


def write_fasta(
    sequences: EncodedSequences,
    path: str | Path,
    names: Sequence[str] | None = None,
) -> None:
    """Write Mimosa sequences while preserving optional FASTA names."""
    labels = (
        tuple(names)
        if names is not None
        else tuple(f"seq{index}" for index in range(len(sequences)))
    )
    if len(labels) != len(sequences):
        raise ValueError("FASTA name count must match sequence count")

    decoder = np.array(["A", "C", "G", "T", "N"], dtype="U1")
    with open(path, "w", encoding="ascii") as handle:
        for name, row in zip(labels, sequences, strict=True):
            clipped = np.clip(row, 0, 4)
            handle.write(f">{name}\n")
            handle.write("".join(decoder[clipped]))
            handle.write("\n")


def write_jstacs_fasta(
    input_fasta: str,
    output_fasta: str,
    position_tag: str = "position",
    value_tag: str = "value",
) -> None:
    """Write FASTA records with Jstacs positional annotations."""
    sequences, names = read_fasta(input_fasta)
    decoder = np.array(["A", "C", "G", "T", "N"], dtype="U1")

    with open(output_fasta, "w", encoding="ascii") as handle:
        for index, (name, row) in enumerate(zip(names, sequences, strict=True)):
            clipped = np.clip(row, 0, 4)
            center = len(row) // 2
            # Jstacs rejects a constant-valued annotation alphabet.
            signal = float(index + 1)
            handle.write(f"> {position_tag}: {center}; {value_tag}: {signal:.1f}\n")
            handle.write("".join(decoder[clipped]))
            handle.write("\n")


def write_meme(
    motifs: Sequence[np.ndarray],
    records: Sequence[tuple[str, int, int]],
    path: str | Path,
) -> None:
    """Write PFMs to a minimal MEME formatted file.

    Each record is ``(name, length, nsites)``; nsites must match every
    PFM column sum, otherwise the reconstructed site counts diverged.
    """
    if len(motifs) != len(records):
        raise ValueError("MEME motif metadata count must match motif count")

    with open(path, "w", encoding="ascii") as handle:
        handle.write("MEME version 4\n\n")
        handle.write("ALPHABET= ACGT\n\n")
        handle.write("strands: + -\n\n")
        handle.write("Background letter frequencies\n")
        handle.write("A 0.25 C 0.25 G 0.25 T 0.25\n\n")

        for pfm, entry in zip(motifs, records, strict=True):
            name, length, nsites = entry
            matrix = np.asarray(pfm, dtype=np.float32)
            if matrix.shape != (4, int(length)):
                raise ValueError(f"PFM for {name} must have shape (4, length)")
            if not np.all(np.isfinite(matrix)):
                raise ValueError(f"PFM for {name} contains non-finite values")
            if not np.allclose(matrix.sum(axis=0), nsites):
                raise ValueError(
                    f"PFM for {name} column sums must all equal nsites={nsites}"
                )
            handle.write(f"MOTIF {name}\n")
            handle.write(
                f"letter-probability matrix: alength= 4 w= {int(length)} "
                f"nsites= {int(nsites)} E= 0\n"
            )
            np.savetxt(handle, matrix.T, fmt=" %.6f")
            handle.write("\n")
