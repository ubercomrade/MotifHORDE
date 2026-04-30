import pytest
import numpy as np
from hordemotifs.functions import (
    pfm_to_pwm, pcm_to_pfm, precision_recall_curve, roc_curve, 
    cut_roc, cut_prc, standardized_pauc, scores_to_frequencies
)
from hordemotifs.ragged import ragged_from_list, RaggedData

@pytest.mark.unit
def test_pfm_to_pwm():
    """Проверка математической корректности pfm_to_pwm (log-odds)."""
    pfm = np.array([[0.8, 0.2], [0.05, 0.2], [0.05, 0.3], [0.1, 0.3]], dtype=np.float32)
    pwm = pfm_to_pwm(pfm)
    # background = 0.25. Для A в поз 0: log((0.8 + 0.0001) / 0.25)
    expected_val = np.log((0.8 + 0.0001) / 0.25)
    assert pwm[0, 0] == pytest.approx(expected_val)

@pytest.mark.unit
def test_pcm_to_pfm():
    """Проверка нормализации pcm_to_pfm."""
    pcm = np.array([[10, 2], [0, 8], [0, 0], [0, 0]], dtype=np.float32)
    pfm = pcm_to_pfm(pcm)
    # column sum of pcm is 10. nuc_pseudo=0.25. (10 + 0.25) / (10 + 1) = 10.25 / 11
    assert pfm[0, 0] == pytest.approx(10.25 / 11.0)
    assert np.allclose(pfm.sum(axis=0), 1.0)

@pytest.mark.unit
def test_precision_recall_curve():
    """Проверка вычисления точности/полноты."""
    y_true = np.array([1, 1, 0, 0], dtype=np.int64)
    y_scores = np.array([0.9, 0.8, 0.7, 0.1], dtype=np.float32)
    prec, rec, thr = precision_recall_curve(y_true, y_scores)
    # При самом высоком пороге (inf) prec=1, rec=0
    assert rec[0] == 0.0
    assert prec[0] == 1.0
    assert rec[-1] == 1.0

@pytest.mark.unit
def test_roc_curve():
    """Проверка вычисления TPR/FPR."""
    y_true = np.array([1, 1, 0, 0], dtype=np.int64)
    y_scores = np.array([0.9, 0.8, 0.7, 0.1], dtype=np.float32)
    tpr, fpr, thr = roc_curve(y_true, y_scores)
    assert tpr[0] == 0.0
    assert fpr[0] == 0.0
    # TP=2, FN=0 -> TPR=1. FP=2, TN=0 -> FPR=1
    assert tpr[-1] == 1.0
    assert fpr[-1] == 1.0

@pytest.mark.unit
def test_curves_edge_cases():
    """Граничные случаи для кривых."""
    # Пустые данные
    prec, rec, thr = precision_recall_curve(np.array([], dtype=np.int64), np.array([], dtype=np.float32))
    assert len(prec) == 1
    # Одинаковые оценки
    y_true = np.array([1, 0], dtype=np.int64)
    y_scores = np.array([0.5, 0.5], dtype=np.float32)
    tpr, fpr, thr = roc_curve(y_true, y_scores)
    assert len(tpr) == 2 # inf + 0.5

@pytest.mark.unit
def test_cut_roc():
    """Обрезка ROC кривой."""
    tpr = np.array([0.0, 0.5, 1.0], dtype=np.float64)
    fpr = np.array([0.0, 0.2, 1.0], dtype=np.float64)
    thr = np.array([10.0, 5.0, 0.0], dtype=np.float64)
    t_c, f_c, s_c = cut_roc(tpr, fpr, thr, score_cutoff=5.0)
    assert s_c[-1] == 5.0
    assert t_c[-1] == 0.5

@pytest.mark.unit
def test_standardized_pauc():
    """Проверка нормализации pAUC к диапазону [0.5, 1]."""
    # Если pauc_raw == pauc_min, результат 0.5
    assert standardized_pauc(0.1, 0.1, 0.2) == 0.5
    # Если pauc_raw == pauc_max, результат 1.0
    assert standardized_pauc(0.2, 0.1, 0.2) == 1.0

@pytest.mark.unit
def test_scores_to_frequencies():
    """Проверка преобразования -log10(FPR)."""
    scores = np.array([10.0, 5.0, 5.0, 0.0], dtype=np.float32)
    ragged = RaggedData(scores, np.array([0, 4], dtype=np.int64))
    freq_ragged = scores_to_frequencies(ragged)
    
    # Всего 4 элемента
    # score=10: rank 1. FPR = 1/4 = 0.25. -log10(0.25) = 0.60206
    # score=5: rank 3. FPR = 3/4 = 0.75. -log10(0.75) = 0.1249
    assert freq_ragged.data[0] == pytest.approx(-np.log10(0.25), abs=1e-5)