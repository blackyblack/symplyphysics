"""
Work of ideal gas in isothermal process
=======================================

The isothermal process of expansion (or compression) of a gas can occur under conditions where heat
exchange between the gas and the external environment is carried out at a constant temperature.
To do this, the heat capacity of the external environment must be large enough, and the expansion
(or compression) process must be slow enough.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.

**Conditions:**

#. The gas is ideal.
#. The temperature of the gas stays constant during the expansion.
"""

from sympy import Eq, solve, log
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import work_is_integral_of_pressure_over_volume as work_law
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation

work = symbols.work
"""
:symbols:`work` done by the ideal gas during the isothermal process.
"""

amount_of_substance = symbols.amount_of_substance
"""
:symbols:`amount_of_substance` of ideal gas.
"""

initial_volume = clone_as_symbol(symbols.volume, subscript="0")
"""
Initial :symbols:`volume` of the gas.
"""

final_volume = clone_as_symbol(symbols.volume, subscript="1")
"""
Final :symbols:`volume` of the gas.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the gas.
"""

law = Eq(
    work,
    amount_of_substance * quantities.molar_gas_constant * temperature *
    log(final_volume / initial_volume),
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive from ideal gas equation

_ideal_gas_eqn = ideal_gas_equation.law.subs({
    ideal_gas_equation.amount_of_substance: amount_of_substance,
    ideal_gas_equation.temperature: temperature,
})

_pressure_expr = solve(_ideal_gas_eqn,
    ideal_gas_equation.pressure)[0].subs(ideal_gas_equation.volume, work_law.volume)

_volume_before = clone_as_symbol(symbols.volume, subscript="0", positive=True)
_volume_after = clone_as_symbol(symbols.volume, subscript="1", positive=True)

_work_expr = work_law.law.rhs.subs({
    work_law.initial_volume: _volume_before,
    work_law.final_volume: _volume_after,
    work_law.pressure(work_law.volume): _pressure_expr
}).doit().simplify()

_work_from_law = law.rhs.subs({
    initial_volume: _volume_before,
    final_volume: _volume_after,
})

assert expr_equals(_work_expr, _work_from_law)


@validate_input(amount_of_substance_=amount_of_substance,
    start_volume_=initial_volume,
    final_volume_=final_volume,
    temperature_=temperature)
@validate_output(work)
def calculate_work(amount_of_substance_: Quantity, start_volume_: Quantity, final_volume_: Quantity,
    temperature_: Quantity) -> Quantity:
    solved = solve(law, work, dict=True)[0][work]
    result_expr = solved.subs({
        amount_of_substance: amount_of_substance_,
        initial_volume: start_volume_,
        final_volume: final_volume_,
        temperature: temperature_
    })
    return Quantity(result_expr)
