"""Input/output helpers for FASTA batches and pipeline motif output."""

from __future__ import annotations

from pathlib import Path
from typing import List, Tuple

import numpy as np
from mimosa.batches import make_sequence_batch, row_values


def read_fasta(path: str | Path):
    """Read a FASTA file and return integer-encoded sequences."""
    trans_table = bytearray([4] * 256)
    for char, code in zip(b"ACGTacgt", [0, 1, 2, 3] * 2, strict=False):
        trans_table[char] = code

    sequences: List[np.ndarray] = []
    with open(path, "r") as handle:
        current_seq_bytes = bytearray()
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_seq_bytes:
                    encoded = np.frombuffer(
                        current_seq_bytes.translate(trans_table), dtype=np.int8
                    ).copy()
                    sequences.append(encoded)
                    current_seq_bytes.clear()
            else:
                current_seq_bytes.extend(line.encode("ascii", errors="ignore"))

        if current_seq_bytes:
            encoded = np.frombuffer(
                current_seq_bytes.translate(trans_table), dtype=np.int8
            ).copy()
            sequences.append(encoded)

    return make_sequence_batch(sequences)


def write_fasta(sequences, path: str | Path) -> None:
    """Write integer-encoded sequences to a FASTA file."""
    decoder = np.array(["A", "C", "G", "T", "N"], dtype="U1")

    if isinstance(sequences, dict) and {"values", "lengths"}.issubset(sequences):
        rows = (
            row_values(sequences, index) for index in range(len(sequences["lengths"]))
        )
    else:
        rows = (np.asarray(row, dtype=np.int8).ravel() for row in sequences)

    with open(path, "w") as handle:
        for index, row in enumerate(rows):
            clipped = np.clip(row, 0, 4)
            handle.write(f">seq{index}\n")
            handle.write("".join(decoder[clipped]))
            handle.write("\n")


def write_jstacs_fasta(
    input_fasta: str,
    output_fasta: str,
    position_tag: str = "position",
    value_tag: str = "value",
) -> None:
    """Write FASTA records with Jstacs positional annotations."""
    sequences = read_fasta(input_fasta)
    decoder = np.array(["A", "C", "G", "T", "N"], dtype="U1")

    with open(output_fasta, "w") as handle:
        for index in range(len(sequences["lengths"])):
            row = row_values(sequences, index)
            clipped = np.clip(row, 0, 4)
            center = int(sequences["lengths"][index]) // 2
            # Jstacs rejects a constant-valued annotation alphabet.
            signal = float(index + 1)
            handle.write(f"> {position_tag}: {center}; {value_tag}: {signal:.1f}\n")
            handle.write("".join(decoder[clipped]))
            handle.write("\n")


def write_meme(
    motifs: List[np.ndarray], info: List[Tuple[str, int]], path: str | Path
) -> None:
    """Write PFMs to a minimal MEME formatted file."""
    with open(path, "w") as handle:
        handle.write("MEME version 4\n\n")
        handle.write("ALPHABET= ACGT\n\n")
        handle.write("strands: + -\n\n")
        handle.write("Background letter frequencies\n")
        handle.write("A 0.25 C 0.25 G 0.25 T 0.25\n\n")

        for pfm, (name, length) in zip(motifs, info, strict=True):
            matrix = np.asarray(pfm, dtype=np.float32)
            if matrix.shape[0] != 4:
                raise ValueError(f"PFM for {name} must have shape (4, length)")
            handle.write(f"MOTIF {name}\n")
            handle.write(
                f"letter-probability matrix: alength= 4 w= {int(length)} nsites= 20 E= 0\n"
            )
            np.savetxt(handle, matrix.T, fmt=" %.6f")
            handle.write("\n")
