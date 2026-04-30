import pytest
import numpy as np
import pandas as pd
from hordemotifs.models import PwmMotif
from hordemotifs.ragged import RaggedData
from tests.test_fixtures import test_pfm, test_pwm, test_ragged_data

@pytest.mark.unit
def test_pwm_motif_initialization(test_pwm, test_pfm):
    """Проверка корректности инициализации PwmMotif."""
    name = "TestMotif"
    length = 6
    motif = PwmMotif(matrix=test_pwm, name=name, length=length, pfm=test_pfm)
    
    assert motif.name == name
    assert motif.length == length
    np.testing.assert_array_almost_equal(motif.matrix, test_pwm)
    np.testing.assert_array_almost_equal(motif.pfm, test_pfm)

@pytest.mark.unit
def test_pwm_motif_properties(test_pwm, test_pfm):
    """Проверка свойств pfm, matrix, name, length."""
    motif = PwmMotif(matrix=test_pwm, name="M", length=6, pfm=test_pfm)
    assert motif.length == 6
    assert motif.name == "M"
    assert motif.matrix.shape == (5, 6)
    assert motif.pfm.shape == (4, 6)

@pytest.mark.unit
def test_scan_forward_strand(test_pwm, test_pfm, test_ragged_data):
    """Сканирование по прямой цепи."""
    motif = PwmMotif(matrix=test_pwm, name="M", length=6, pfm=test_pfm)
    results = motif.scan(test_ragged_data, strand="+")
    
    assert isinstance(results, RaggedData)
    assert results.num_sequences == test_ragged_data.num_sequences
    for i in range(results.num_sequences):
        expected_len = test_ragged_data.get_length(i) - motif.length + 1
        assert results.get_length(i) == expected_len

@pytest.mark.unit
def test_scan_reverse_strand(test_pwm, test_pfm, test_ragged_data):
    """Сканирование по обратной цепи."""
    motif = PwmMotif(matrix=test_pwm, name="M", length=6, pfm=test_pfm)
    results = motif.scan(test_ragged_data, strand="-")
    assert results.num_sequences == test_ragged_data.num_sequences

@pytest.mark.unit
def test_scan_both_strands(test_pwm, test_pfm, test_ragged_data):
    """Сканирование обеих цепей."""
    motif = PwmMotif(matrix=test_pwm, name="M", length=6, pfm=test_pfm)
    results = motif.scan(test_ragged_data, strand="both")
    # Должно быть 2 * N последовательностей
    assert results.num_sequences == 2 * test_ragged_data.num_sequences

@pytest.mark.unit
def test_scan_best_strand(test_pwm, test_pfm, test_ragged_data):
    """Сканирование лучшей цепи (max per position)."""
    motif = PwmMotif(matrix=test_pwm, name="M", length=6, pfm=test_pfm)
    results = motif.scan(test_ragged_data, strand="best")
    assert results.num_sequences == test_ragged_data.num_sequences

@pytest.mark.unit
def test_best_scores(test_pwm, test_pfm, test_ragged_data):
    """Лучшие оценки для последовательностей."""
    motif = PwmMotif(matrix=test_pwm, name="M", length=6, pfm=test_pfm)
    scores = motif.best_scores(test_ragged_data, strand="best")
    assert len(scores) == test_ragged_data.num_sequences
    assert scores.dtype == np.float32

@pytest.mark.unit
def test_threshold_table_computation(test_pwm, test_pfm, test_ragged_data):
    """Вычисление таблицы порогов."""
    motif = PwmMotif(matrix=test_pwm, name="M", length=6, pfm=test_pfm)
    table = motif.get_threshold_table(test_ragged_data)
    assert table.shape[1] == 2
    assert np.all(np.diff(table[:, 0]) <= 0) # Сортировка по убыванию score

@pytest.mark.unit
def test_score_to_frequency(test_pwm, test_pfm, test_ragged_data):
    """Преобразование score в частоту."""
    motif = PwmMotif(matrix=test_pwm, name="M", length=6, pfm=test_pfm)
    motif.get_threshold_table(test_ragged_data)
    max_score = motif.threshold_table[0, 0]
    freq = motif._score_to_frequency(max_score)
    assert freq >= 0

@pytest.mark.unit
def test_frequency_to_score(test_pwm, test_pfm, test_ragged_data):
    """Преобразование частоты (FPR) в score."""
    motif = PwmMotif(matrix=test_pwm, name="M", length=6, pfm=test_pfm)
    motif.get_threshold_table(test_ragged_data)
    # 0.01 -> -log10(0.01) = 2
    score = motif._frequency_to_score(0.01)
    assert isinstance(score, float)

@pytest.mark.unit
def test_get_sites_best_mode(test_pwm, test_pfm, test_ragged_data):
    """Поиск лучшего сайта."""
    motif = PwmMotif(matrix=test_pwm, name="M", length=6, pfm=test_pfm)
    df = motif.get_sites(test_ragged_data, mode="best")
    assert isinstance(df, pd.DataFrame)
    assert len(df) <= test_ragged_data.num_sequences
    assert "score" in df.columns

@pytest.mark.unit
def test_get_sites_threshold_mode(test_pwm, test_pfm, test_ragged_data):
    """Поиск сайтов по порогу."""
    motif = PwmMotif(matrix=test_pwm, name="M", length=6, pfm=test_pfm)
    motif.get_threshold_table(test_ragged_data)
    df = motif.get_sites(test_ragged_data, mode="threshold", fpr_threshold=0.1)
    assert isinstance(df, pd.DataFrame)