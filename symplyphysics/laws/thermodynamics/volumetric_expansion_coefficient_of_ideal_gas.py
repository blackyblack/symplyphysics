from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import volumetric_coefficient_of_thermal_expansion as coef_def
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation

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

# Derive from ideal gas equation and definition of volumetric expansion coefficient

_volume_expr = solve(ideal_gas_equation.law,
    ideal_gas_equation.volume)[0].subs(ideal_gas_equation.temperature, temperature)

_coef_expr = coef_def.definition.rhs.subs(coef_def.temperature,
    temperature).subs(coef_def.volume(temperature), _volume_expr).doit()

assert expr_equals(_coef_expr, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(temperature_=temperature)
@validate_output(volumetric_expansion_coefficient)
def calculate_volumetric_expansion_coefficient(temperature_: Quantity) -> Quantity:
    result = law.rhs.subs(temperature, temperature_)
    return Quantity(result)
