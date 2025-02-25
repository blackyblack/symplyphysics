from __future__ import annotations


class Probability(float):

    def __new__(cls, value: float) -> Probability:
        if value < 0 or value > 1.0:
            raise AttributeError(f"Probability value should be in range [0..1], got '{value}'")
        return float.__new__(cls, value)
