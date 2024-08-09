"""
Volumetric expansion coefficient of ideal gas
=============================================

The isobaric volumetric expansion coefficient of an ideal gas is the inverse of its temperature.

**Conditions:**

#. Pressure of the gas remains constant during the expansion.
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import volumetric_coefficient_of_thermal_expansion as coef_def
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation

volumetric_expansion_coefficient = Symbol("volumetric_expansion_coefficient", 1 / units.temperature)
r"""
:doc:`Volumetric expansion coefficient <definitions.volumetric_coefficient_of_thermal_expansion>` of the gas.

Symbol:
    :code:`alpha_V`

Latex:
    :math:`\alpha_V`
"""

temperature = symbols.thermodynamics.temperature
"""
:attr:`~symplyphysics.symbols.thermodynamics.temperature` of the gas.

Symbol:
    :code:`T`
"""

law = Eq(volumetric_expansion_coefficient, 1 / temperature)
r"""
:code:`alpha_V = 1 / T`

Latex:
    .. math::
        \alpha_V = \frac{1}{T}
"""

# Derive from ideal gas equation and definition of volumetric expansion coefficient

_volume_expr = solve(ideal_gas_equation.law,
    ideal_gas_equation.volume)[0].subs(ideal_gas_equation.temperature, temperature)

_coef_expr = coef_def.definition.rhs.subs(coef_def.temperature,
    temperature).subs(coef_def.volume(temperature), _volume_expr).doit()

assert expr_equals(_coef_expr, law.rhs)


@validate_input(temperature_=temperature)
@validate_output(volumetric_expansion_coefficient)
def calculate_volumetric_expansion_coefficient(temperature_: Quantity) -> Quantity:
    result = law.rhs.subs(temperature, temperature_)
    return Quantity(result)
