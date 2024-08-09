"""
Gas pressure change from temperature
====================================

The change in pressure of ideal gas depends on final_temperature, standard pressure and thermal coefficient.

**Conditions:**

#. Gas is ideal.
#. The heating of the gas is isochoric.
"""

from sympy import Eq, solve
from symplyphysics import (
    symbols,
    clone_symbol,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    quantities,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import pressure_and_temperature_in_isochoric_process as isochoric_law

pressure_change = Symbol("pressure", units.pressure)
r"""
Finite change in pressure inside the gas.

Symbol:
    :code:`dp`

Latex:
    :math:`\Delta p`
"""

initial_pressure = Symbol("initial_pressure", units.pressure)
r"""
Pressure at :attr:`standard conditions <symplyphysics.quantities.standard_conditions_temperature>`.

Symbol:
    :code:`p_0`

Latex:
    :math:`p_0`
"""

thermal_coefficient = 1 / quantities.standard_conditions_temperature
r"""
Isochoric thermal pressure coefficient is equal to the inversed temperature of standard conditions.

Symbol:
    :code:`beta_V`

Latex:
    :math:`\beta_V`
"""

final_temperature = clone_symbol(symbols.thermodynamics.temperature, "final_temperature")
"""
:attr:`~symplyphysics.symbols.thermodynamics.temperature` of the gas.

Symbol:
    :code:`T`
"""

law = Eq(pressure_change, initial_pressure * (thermal_coefficient * final_temperature - 1))
r"""
:code:`dp = p_0 * (beta_V * T - 1)`

Latex:
    .. math::
        \Delta p = p_0 \left( \beta_V T - 1 \right)
"""

# Derive law from Charles' law of isochoric heating of gas.

_isochoric_eqn = isochoric_law.law.subs({
    isochoric_law.initial_pressure: initial_pressure,
    isochoric_law.final_pressure: initial_pressure + pressure_change,
    isochoric_law.initial_temperature: quantities.standard_conditions_temperature,
    isochoric_law.final_temperature: final_temperature,
})

_pressure_change_expr = solve(_isochoric_eqn, pressure_change)[0]

assert expr_equals(_pressure_change_expr, law.rhs)


@validate_input(standard_pressure_=initial_pressure, temperature_=final_temperature)
@validate_output(pressure_change)
def calculate_pressure_change(standard_pressure_: Quantity, temperature_: Quantity) -> Quantity:
    result_expr = law.rhs.subs({initial_pressure: standard_pressure_, final_temperature: temperature_})
    return Quantity(result_expr)
