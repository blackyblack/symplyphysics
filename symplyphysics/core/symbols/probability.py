from dataclasses import dataclass


@dataclass(frozen=True)
class Probability:
    value: float

    def __post_init__(self):
        if self.value < 0 or self.value > 1.0:
            raise AttributeError(f"Probability value should be in range [0..1], got '{self.value}'")
