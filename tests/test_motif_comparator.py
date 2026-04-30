import pytest
import numpy as np
import pandas as pd
from hordemotifs.comparison import MotifComparator
from hordemotifs.models import PwmMotif
from tests.test_fixtures import test_pfm, test_pwm, test_ragged_data

@pytest.fixture
def test_motifs(test_pwm, test_pfm):
    m1 = PwmMotif(matrix=test_pwm, name="M1", length=6, pfm=test_pfm)
    m2 = PwmMotif(matrix=test_pwm, name="M2", length=6, pfm=test_pfm)
    return [m1], [m2]

@pytest.mark.unit
def test_jaccard_metric_calculation(test_motifs, test_ragged_data):
    """Вычисление CJ (Jaccard-like) индекса."""
    m1_list, m2_list = test_motifs
    comparator = MotifComparator(metric='cj', n_permutations=0)
    
    freq1 = m1_list[0].get_frequencies(test_ragged_data, strand="+")
    freq2 = m2_list[0].get_frequencies(test_ragged_data, strand="+")
    
    score, offset = comparator._compute_metric(freq1, freq2, search_range=5, metric='cj')
    assert score > 0.99
    assert offset == 0

@pytest.mark.unit
def test_overlap_metric_calculation(test_motifs, test_ragged_data):
    """Вычисление CO (Overlap) коэффициента."""
    m1_list, m2_list = test_motifs
    comparator = MotifComparator(metric='co', n_permutations=0)
    
    freq1 = m1_list[0].get_frequencies(test_ragged_data, strand="+")
    freq2 = m2_list[0].get_frequencies(test_ragged_data, strand="+")
    
    score, offset = comparator._compute_metric(freq1, freq2, search_range=5, metric='co')
    assert score > 0.99
    assert offset == 0

@pytest.mark.unit
def test_single_compare_no_stats(test_motifs, test_ragged_data):
    """Сравнение без статистики."""
    m1_list, m2_list = test_motifs
    comparator = MotifComparator(n_permutations=0)
    result = comparator._single_compare(m1_list[0], m2_list[0], test_ragged_data, stat=False)
    
    assert "score" in result
    assert "offset" in result
    assert "p-value" not in result

@pytest.mark.unit
def test_single_compare_with_stats(test_motifs, test_ragged_data):
    """Сравнение с перестановочными тестами (детерминировано)."""
    m1_list, m2_list = test_motifs
    # Уменьшим число перестановок для быстроты теста
    comparator = MotifComparator(n_permutations=10, seed=42, n_jobs=1)
    result = comparator._single_compare(m1_list[0], m2_list[0], test_ragged_data, stat=True)
    
    assert "p-value" in result
    assert "z-score" in result
    assert 0 <= result["p-value"] <= 1

@pytest.mark.unit
def test_compare_multiple_motifs(test_motifs, test_ragged_data):
    """Сравнение нескольких мотивов."""
    m1_list, m2_list = test_motifs
    comparator = MotifComparator(n_permutations=0)
    df = comparator.compare(m1_list, m2_list, sequences=test_ragged_data)
    
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 1
    assert df.loc[0, "query"] == "M1"
    assert df.loc[0, "target"] == "M2"

@pytest.mark.unit
def test_strand_selection(test_motifs, test_ragged_data):
    """Выбор лучшей ориентации (++, +-)."""
    m1_list, m2_list = test_motifs
    comparator = MotifComparator(n_permutations=0)
    result = comparator._single_compare(m1_list[0], m2_list[0], test_ragged_data, stat=False)
    
    assert result["strand_pair"] in ["++", "+-"]

@pytest.mark.unit
def test_filtering_by_score(test_motifs, test_ragged_data):
    """Фильтрация по score."""
    m1_list, m2_list = test_motifs
    comparator = MotifComparator(n_permutations=0, filter_type='score', filter_threshold=0.9)
    df = comparator.compare(m1_list, m2_list, sequences=test_ragged_data)
    assert len(df) == 1
    
    comparator = MotifComparator(n_permutations=0, filter_type='score', filter_threshold=1.1)
    df = comparator.compare(m1_list, m2_list, sequences=test_ragged_data)
    assert len(df) == 0

@pytest.mark.unit
def test_filtering_by_pvalue(test_motifs, test_ragged_data):
    """Фильтрация по p-value."""
    m1_list, m2_list = test_motifs
    comparator = MotifComparator(n_permutations=10, seed=42, filter_type='p-value', filter_threshold=1.0, n_jobs=1)
    df = comparator.compare(m1_list, m2_list, sequences=test_ragged_data)
    assert len(df) == 1