class Probability(object):
    def __init__(self, value: float = 0):
        self.value = value

    @property
    def value(self) -> float:
        return self._value

    @value.setter
    def value(self, value: float):
        if value < 0 or value > 1.0:
            raise AttributeError(f"Probability value should be in range [0..1], got '{value}'")
        self._value = value
