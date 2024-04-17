from sympy import Eq, symbols, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import (
    linear_coefficient_of_thermal_expansion as linear_def,
    volumetric_coefficient_of_thermal_expansion as volumetric_def,
)

# Description
## Coefficients of thermal expansion describe how the size of an object changes with a change in temperature
## at a constant pressure. In isotropic materials, the volumetric coefficient is three times the linear one.

# Law: beta = 3 * alpha
## beta - volumetric coefficient of thermal expansion
## alpha - linear coefficient of thermal expansion

volumetric_expansion_coefficient = Symbol("volumetric_expansion_coefficient", 1 / units.temperature)
linear_expansion_coefficient = Symbol("linear_expansion_coefficient", 1 / units.temperature)

law = Eq(volumetric_expansion_coefficient, 3 * linear_expansion_coefficient)

# Derive from their definitions for a cube-sized object

_volume_change = symbols("volume_change")
_initial_length, _length_change = symbols("initial_length length_change")
_temperature_change = symbols("temperature_change")

_final_length_expr = _initial_length + _length_change
_initial_volume_expr = _initial_length**3
_final_volume_expr = _final_length_expr**3

_volumetric_expr = volumetric_def.definition.rhs.subs({
    volumetric_def.volume(volumetric_def.temperature).diff(volumetric_def.temperature): _volume_change / _temperature_change,
    _volume_change: _final_volume_expr - _initial_volume_expr,
    volumetric_def.volume(volumetric_def.temperature): _initial_volume_expr,
})

# Note: if dimension changes are not small, the third and fourth terms of the following series expansion
# can be taken into account.
_volumetric_expr_series = _volumetric_expr.series(_length_change, 0, 2).removeO()

_linear_eqn = linear_def.definition.subs({
    linear_def.length(linear_def.temperature).diff(linear_def.temperature): _length_change / _temperature_change,
    linear_def.length(linear_def.temperature): _initial_length,
    linear_def.linear_expansion_coefficient: linear_expansion_coefficient,
})

_volumetric_expr_derived = solve(
    (Eq(volumetric_expansion_coefficient, _volumetric_expr_series), _linear_eqn),
    (volumetric_expansion_coefficient, _length_change),
    dict=True,
)[0][volumetric_expansion_coefficient]

assert expr_equals(_volumetric_expr_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(linear_expansion_coefficient_=linear_expansion_coefficient)
@validate_output(volumetric_expansion_coefficient)
def calculate_volumetric_expansion_coefficient(linear_expansion_coefficient_: Quantity) -> Quantity:
    result = law.rhs.subs(linear_expansion_coefficient, linear_expansion_coefficient_)
    return Quantity(result)
