from typing import Self


class Probability(float):

    def __new__(cls, value: float) -> Self:
        if value < 0 or value > 1.0:
            raise AttributeError(f"Probability value should be in range [0..1], got '{value}'")
        return float.__new__(cls, value)
