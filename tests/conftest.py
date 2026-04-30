import pytest
import numpy as np
import os
import tempfile
from hordemotifs.ragged import ragged_from_list

@pytest.fixture
def tmp_dir():
    """Создает временную директорию для тестов."""
    with tempfile.TemporaryDirectory() as tmp:
        yield tmp

@pytest.fixture
def sample_fasta(tmp_dir):
    """Создает тестовый FASTA файл."""
    path = os.path.join(tmp_dir, "test.fa")
    with open(path, "w") as f:
        f.write(">seq0\nACGTACGTACGTACGTACGT\n")
        f.write(">seq1\nGATTACAGATTACAGATTA\n")
    return path

@pytest.fixture
def sample_meme(tmp_dir):
    """Создает тестовый MEME файл."""
    path = os.path.join(tmp_dir, "test.meme")
    content = """MEME version 4
ALPHABET= ACGT
strands: + -

MOTIF TestMotif1
letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0
  0.800000  0.050000  0.050000  0.100000
  0.100000  0.800000  0.050000  0.050000
  0.100000  0.050000  0.800000  0.050000
  0.100000  0.100000  0.100000  0.700000
  0.200000  0.200000  0.300000  0.300000
  0.250000  0.250000  0.250000  0.250000
"""
    with open(path, "w") as f:
        f.write(content)
    return path

@pytest.fixture(autouse=True)
def set_seed():
    """Устанавливает seed для детерминизма."""
    np.random.seed(42)