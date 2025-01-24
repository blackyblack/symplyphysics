"""
Gas pressure change from temperature
====================================

The change in pressure of ideal gas depends on temperature, standard pressure and thermal coefficient.

**Conditions:**

#. Gas is ideal.
#. The heating of the gas is isochoric.

..
    TODO find link
"""

from sympy import Eq, solve
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
    quantities,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import pressure_and_temperature_in_isochoric_process as isochoric_law

pressure_change = clone_as_symbol(symbols.pressure, display_symbol="Delta(p)", display_latex="\\Delta p")
"""
Finite change in :symbols:`pressure` inside the gas.
"""

initial_pressure = clone_as_symbol(symbols.pressure, subscript="0")
"""
:symbols:`pressure` at :attr:`standard conditions <symplyphysics.quantities.standard_conditions_temperature>`.
"""

thermal_coefficient = 1 / quantities.standard_conditions_temperature
"""
Isochoric thermal pressure coefficient equal to the inverse of
:attr:`standard conditions temperature <symplyphysics.quantities.standard_conditions_temperature>`.
"""

final_temperature = symbols.temperature
"""
:symbols:`temperature` of the gas.
"""

law = Eq(pressure_change, initial_pressure * (thermal_coefficient * final_temperature - 1))
"""
:laws:symbol::

:laws:latex::
"""

# Derive law from Guy-Lussac's law of isochoric heating of gas.

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
    result_expr = law.rhs.subs({
        initial_pressure: standard_pressure_,
        final_temperature: temperature_
    })
    return Quantity(result_expr)
