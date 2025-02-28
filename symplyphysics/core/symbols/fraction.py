from __future__ import annotations


class Fraction(float):

    def __new__(cls, value: float) -> Fraction:
        if value < 0 or value > 1.0:
            raise AttributeError(f"Fraction value should be in range [0..1], got '{value}'")
        return float.__new__(cls, value)
