import pytest
import os
import numpy as np
import subprocess
from unittest.mock import MagicMock, patch
from hordemotifs.discovery import StremeDiscoveryTool
from hordemotifs.models import PwmMotif

@pytest.mark.integration
def test_streme_discovery_initialization():
    """Проверка инициализации StremeDiscoveryTool."""
    tool = StremeDiscoveryTool(nmotifs=3)
    assert tool.name == "streme"
    assert tool.nmotifs == 3

@pytest.mark.integration
@patch("subprocess.run")
def test_streme_discovery_with_mock(mock_run, tmp_dir, sample_fasta, sample_meme):
    """Использование mock для subprocess.run и проверка парсинга."""
    # Подготавливаем mock для вызова streme
    with open(sample_meme, "rb") as f:
        mock_stdout = f.read()
    
    mock_run.return_value = MagicMock(stdout=mock_stdout, returncode=0)
    
    tool = StremeDiscoveryTool(nmotifs=1)
    # Вызываем discover
    motifs = tool.discover(
        foreground=sample_fasta,
        background=sample_fasta,
        output_dir=tmp_dir,
        number_of_motifs=1,
        length=6
    )
    
    # Проверяем вызов
    assert mock_run.called
    args = mock_run.call_args[0][0]
    assert "streme" in args
    assert "--w" in args
    
    # Проверяем результат
    assert len(motifs) > 0
    assert isinstance(motifs[0], PwmMotif)
    assert motifs[0].name == "Streme-1"
    assert motifs[0].length == 6
    assert motifs[0].pfm.shape == (4, 6)

@pytest.mark.unit
def test_pwm_creation_from_pfm():
    """Проверка создания PWM из PFM в контексте discovery."""
    from hordemotifs.functions import pfm_to_pwm
    pfm = np.array([[0.8, 0.2], [0.05, 0.2], [0.05, 0.3], [0.1, 0.3]], dtype=np.float32)
    pwm = pfm_to_pwm(pfm)
    pwm_ext = np.concatenate((pwm, np.min(pwm, axis=0, keepdims=True)), axis=0)
    
    assert pwm_ext.shape == (5, 2)
    assert np.all(pwm_ext[4, :] == np.min(pwm, axis=0))

@pytest.mark.integration
@patch("subprocess.run")
def test_streme_no_motifs_found(mock_run, tmp_dir, sample_fasta):
    """Случай, когда STREME не нашел мотивов (пустой вывод)."""
    mock_run.return_value = MagicMock(stdout=b"", returncode=0)
    
    tool = StremeDiscoveryTool(nmotifs=1)
    motifs = tool.discover(
        foreground=sample_fasta,
        background=sample_fasta,
        output_dir=tmp_dir,
        number_of_motifs=1,
        length=6
    )
    assert len(motifs) == 0