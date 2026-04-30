import numpy as np
import pytest
from hordemotifs.ragged import ragged_from_list
from hordemotifs.functions import pfm_to_pwm

def get_test_pfm():
    """
    Возвращает детерминированную PFM (4, 6).
    """
    pfm = np.array([
        [0.8, 0.1, 0.1, 0.1, 0.2, 0.25], # A
        [0.05, 0.8, 0.05, 0.1, 0.2, 0.25], # C
        [0.05, 0.05, 0.8, 0.1, 0.3, 0.25], # G
        [0.1, 0.05, 0.05, 0.7, 0.3, 0.25]  # T
    ], dtype=np.float32)
    return pfm

def get_test_pwm():
    """
    Вычисляет PWM из PFM с расширением (min score).
    """
    pfm = get_test_pfm()
    pwm = pfm_to_pwm(pfm)
    # Добавляем 5-ю строку с минимумами (как ожидает PwmMotif)
    pwm_ext = np.concatenate([pwm, np.min(pwm, axis=0, keepdims=True)], axis=0)
    return pwm_ext

def get_test_sequences():
    """
    Список из 5 коротких последовательностей.
    A=0, C=1, G=2, T=3, N=4
    """
    seq_strings = [
        "ACGTACGTACGTACGTACGT", # 0: повтор ACGT
        "GATTACAGATTACAGATTA",  # 1: GATTACA
        "TTTTTTTTTTTTTTTTTTTT", # 2: Poly-T
        "AAAAAACCCCCCCGGGGGG",  # 3: Structured
        "ATGCTAGCTAGCTAGCTAGC"  # 4: Random-ish
    ]
    
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    seqs = []
    for s in seq_strings:
        seqs.append(np.array([mapping[c] for c in s], dtype=np.int8))
    return seqs

def get_test_ragged_data():
    """
    RaggedData из тестовых последовательностей.
    """
    return ragged_from_list(get_test_sequences())

@pytest.fixture
def test_pfm():
    return get_test_pfm()

@pytest.fixture
def test_pwm():
    return get_test_pwm()

@pytest.fixture
def test_sequences():
    return get_test_sequences()

@pytest.fixture
def test_ragged_data():
    return get_test_ragged_data()