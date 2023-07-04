from sympy.physics import units
from sympy.physics.units import convert_to
from .quantities import Quantity


class Celsius:
    CELSIUS_TO_KELVIN_OFFSET = 273.15
    value: float

    def __init__(self, value: float = 0):
        self.value = value

    def __str__(self) -> str:
        return str(self.value)

    @property
    def dimension_str(self) -> str:
        return "celsius"


def to_kelvin(value: Celsius) -> float:
    return value.value + Celsius.CELSIUS_TO_KELVIN_OFFSET


def to_kelvin_quantity(value: Celsius) -> Quantity:
    return Quantity(to_kelvin(value) * units.kelvin)


# Here we allow negative Kelvin temperatures, but it does not matter for us. It's
# up to the user to verify that Kelvin temperature is not below zero.
def from_kelvin(value: float) -> Celsius:
    return Celsius(value - Celsius.CELSIUS_TO_KELVIN_OFFSET)


def from_kelvin_quantity(value: Quantity) -> Celsius:
    kelvin_value = float(convert_to(value, units.kelvin).subs(units.kelvin, 1).evalf())
    return from_kelvin(kelvin_value)
