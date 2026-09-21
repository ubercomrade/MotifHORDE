from __future__ import annotations

import numpy as np
import pytest
from mimosa import EncodedSequences, PWM


@pytest.fixture
def sequences() -> EncodedSequences:
    return EncodedSequences.from_strings(
        ["ACGTACGTACGT", "AAAAACCCCCGG", "GGGGTTTTAAAA", "TACGTACGTACG"]
    )


@pytest.fixture
def pwm_model() -> PWM:
    weights = np.array(
        [
            [3.0, -2.0, -2.0, 3.0],
            [-2.0, 3.0, -2.0, -2.0],
            [-2.0, -2.0, 3.0, -2.0],
            [-2.0, -2.0, -2.0, 3.0],
            [-3.0, -3.0, -3.0, -3.0],
        ],
        dtype=np.float32,
    )
    return PWM("Fixture", weights, (0.25, 0.25, 0.25, 0.25))
