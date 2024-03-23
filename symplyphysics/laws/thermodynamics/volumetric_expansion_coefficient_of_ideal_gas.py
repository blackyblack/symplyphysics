from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols,
)

# Description
## The isobaric volumetric expansion coefficient of an ideal gas is the inverse of its temperature.

# Law: alpha_V = 1 / T
## alpha_V = volumetric expansion coefficient
## T - temperature

# Conditions
## - Pressure remains constant during expansion.

volumetric_expansion_coefficient = Symbol("volumetric_expansion_coefficient", 1 / units.temperature)
temperature = symbols.thermodynamics.temperature

law = Eq(volumetric_expansion_coefficient, 1 / temperature)

# TODO: derive from ideal gas equation and definition of volumetric expansion coefficient


def print_law() -> str:
    return print_expression(law)


@validate_input(temperature_=temperature)
@validate_output(volumetric_expansion_coefficient)
def calculate_volumetric_expansion_coefficient(temperature_: Quantity) -> Quantity:
    result = law.rhs.subs(temperature, temperature_)
    return Quantity(result)
