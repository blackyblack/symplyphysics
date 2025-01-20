"""
Volumetric expansion coefficient of ideal gas
=============================================

The isobaric volumetric expansion coefficient of an ideal gas is the inverse of its temperature.

**Conditions:**

#. Pressure of the gas remains constant during the expansion.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Thermal_expansion#In_gases>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import volumetric_coefficient_of_thermal_expansion as coef_def
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation

volumetric_expansion_coefficient = clone_as_symbol(symbols.thermal_expansion_coefficient,
    subscript="V")
"""
Volumetric :symbols:`thermal_expansion_coefficient` of the material. Also see
:doc:`Volumetric expansion coefficient <definitions.volumetric_coefficient_of_thermal_expansion>`.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the gas.
"""

law = Eq(volumetric_expansion_coefficient, 1 / temperature)
"""
:laws:symbol::

:laws:latex::
"""

# Derive from ideal gas equation and definition of volumetric expansion coefficient

_volume_expr = solve(ideal_gas_equation.law,
    ideal_gas_equation.volume)[0].subs(ideal_gas_equation.temperature, temperature)

_coef_expr = coef_def.definition.rhs.subs(coef_def.temperature, temperature)
_coef_expr = _coef_expr.subs(coef_def.volume(temperature, coef_def.parameters), _volume_expr).doit()

assert expr_equals(_coef_expr, law.rhs)


@validate_input(temperature_=temperature)
@validate_output(volumetric_expansion_coefficient)
def calculate_volumetric_expansion_coefficient(temperature_: Quantity) -> Quantity:
    result = law.rhs.subs(temperature, temperature_)
    return Quantity(result)
