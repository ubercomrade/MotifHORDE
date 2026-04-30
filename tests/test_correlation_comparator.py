import pytest
import numpy as np
import pandas as pd
from hordemotifs.comparison import CorrelationComparator
from hordemotifs.models import PwmMotif
from hordemotifs.ragged import RaggedData, ragged_from_list


def test_correlation_comparator():
    """Test the CorrelationComparator implementation."""
    
    # Create dummy sequences as RaggedData
    sequences_data = [
        np.array([0, 1, 2, 3, 0, 1, 2, 3], dtype=np.int8),  # A, C, G, T, A, C, G, T
        np.array([1, 2, 3, 0, 1, 2], dtype=np.int8),        # C, G, T, A, C, G
        np.array([2, 3, 0, 1, 2, 3, 0], dtype=np.int8),    # G, T, A, C, G, T, A
    ]
    sequences = ragged_from_list(sequences_data)
    
    # Create two dummy MotifModel objects (PwmMotif instances)
    # Create simple PFM matrices
    pfm1 = np.array([
        [0.25, 0.25, 0.25, 0.25],  # Equal probabilities for all nucleotides
        [0.5,  0.25, 0.1,  0.15],  # More A
        [0.1,  0.1,  0.7,  0.1],   # More G
        [0.2,  0.3,  0.3,  0.2],   # Balanced
        [0.8,  0.1,  0.05, 0.05],  # Mostly A
    ], dtype=np.float32).T
    
    pfm2 = np.array([
        [0.3,  0.2,  0.3,  0.2],   # Balanced
        [0.1,  0.7,  0.1,  0.1],   # More C
        [0.25, 0.25, 0.25, 0.25],  # Equal probabilities
        [0.15, 0.1,  0.25, 0.5],   # More T
        [0.05, 0.05, 0.1,  0.8],   # Mostly T
    ], dtype=np.float32).T
    
    # Convert PFM to PWM (log odds)
    from hordemotifs.functions import pfm_to_pwm
    pwm1_extended = np.vstack([pfm_to_pwm(pfm1), np.min(pfm_to_pwm(pfm1), axis=0)])
    pwm2_extended = np.vstack([pfm_to_pwm(pfm2), np.min(pfm_to_pwm(pfm2), axis=0)])
    
    motif1 = PwmMotif(
        matrix=pwm1_extended,
        name="motif1",
        length=5,
        pfm=pfm1
    )
    
    motif2 = PwmMotif(
        matrix=pwm2_extended,
        name="motif2",
        length=5,
        pfm=pfm2
    )
    
    # Instantiate CorrelationComparator
    comparator = CorrelationComparator(name="test_correlation")
    
    # Call the compare method
    result = comparator.compare([motif1], [motif2], sequences=sequences)
    
    # Assertions
    assert isinstance(result, pd.DataFrame), "Result should be a pandas DataFrame"
    
    # Check that the result has the expected columns
    expected_columns = ['query', 'target', 'motif1_strand', 'motif2_strand', 'strand_pair', 'score', 'offset', 'search_range']
    for col in expected_columns:
        assert col in result.columns, f"Column '{col}' should be present in the result"
    
    # Check that 'score' column contains values between -1 and 1 (Pearson correlation)
    assert all(-1.0 <= score <= 1.0 for score in result['score']), \
           "All scores should be between -1 and 1 (Pearson correlation range)"
    
    # Check that 'offset' values are integers within reasonable range
    assert all(isinstance(offset, (int, np.integer)) for offset in result['offset']), \
           "All offsets should be integers"
    
    # Check that 'strand_pair' values are correctly reported
    expected_strand_pairs = ['++', '+-', '-+', '--']
    assert all(strand_pair in ['++', '+-'] for strand_pair in result['strand_pair']), \
           "All strand pairs should be either '++' or '+-' (only forward and reverse complement)"
    
    # Additional check: ensure we have one result row for the one vs one motif comparison
    assert len(result) == 1, f"Expected 1 result row, got {len(result)}"
    
    print("All tests passed!")
    print(f"Result DataFrame:\n{result}")
    print(f"Score: {result.iloc[0]['score']}")
    print(f"Offset: {result.iloc[0]['offset']}")
    print(f"Strand pair: {result.iloc[0]['strand_pair']}")