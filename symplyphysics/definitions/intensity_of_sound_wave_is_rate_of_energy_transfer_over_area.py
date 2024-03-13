from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The intensity of a sound wave at a surface is the average rate per unit area at which
## energy is transferred by the wave through or onto the surface.

# Law: I = P / A
## I - intensity of sound wave
## P - rate of energy transfer (power) of sound wave
## A - surface area

intensity = Symbol("intensity", units.power / units.area)
power = Symbol("power", units.power)
area = Symbol("area", units.area)

definition = Eq(intensity, power / area)


def print_law() -> str:
    return print_expression(definition)


@validate_input(power_=power, area_=area)
@validate_output(intensity)
def calculate_intensity(power_: Quantity, area_: Quantity) -> Quantity:
    result = definition.rhs.subs({power: power_, area: area_})
    return Quantity(result)
