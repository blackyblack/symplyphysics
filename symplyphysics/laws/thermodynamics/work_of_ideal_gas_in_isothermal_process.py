r"""
Work of ideal gas in isothermal process
=======================================

The isothermal process of expansion (or compression) of a gas can occur under conditions where heat
exchange between the gas and the external environment is carried out at a constant temperature.
To do this, the heat capacity of the external environment must be large enough, and the expansion
(or compression) process must be slow enough.

**Notation:**

#. :math:`R` is the molar gas constant.
"""

from sympy import Eq, solve, log, Symbol as SymSymbol
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import work_is_volume_integral_of_pressure as work_law
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation

work = Symbol("work", units.energy)
"""
Work done by the ideal gas during the isothermal process.

Symbol:
    :code:`W`
"""

amount_of_substance = Symbol("amount_of_substance", units.amount_of_substance)
"""
Amount of ideal gas.

Symbol:
    :code:`n`
"""

initial_volume = Symbol("initial_volume", units.volume)
r"""
Initial volume of the gas.

Symbol:
    :code:`V_0`

Latex:
    :math:`V_0`
"""

final_volume = Symbol("final_volume", units.volume)
"""
Final volume of the gas.

Symbol:
    :code:`V_1`

Latex:
    :math:`V_1`
"""

temperature = symbols.thermodynamics.temperature
"""
:attr:`~symplyphysics.symbols.thermodynamics.temperature` of the gas, constant throughout the process.

Symbol:
    :code:`T`
"""

law = Eq(work, amount_of_substance * units.molar_gas_constant * temperature *
    log(final_volume / initial_volume))
r"""
:code:`W = n * R * T * log(V_1 / V_0)`

Latex:
    .. math::
        W = n R T \log \frac{V_1}{V_0}
"""

# Derive from ideal gas equation

_ideal_gas_eqn = ideal_gas_equation.law.subs({
    ideal_gas_equation.mole_count: amount_of_substance,
    ideal_gas_equation.temperature: temperature,
})

_pressure_expr = solve(_ideal_gas_eqn,
    ideal_gas_equation.pressure)[0].subs(ideal_gas_equation.volume, work_law.volume)

_volume_before = SymSymbol("volume_before", positive=True)
_volume_after = SymSymbol("volume_after", positive=True)

_work_expr = work_law.law.rhs.subs({
    work_law.volume_before: _volume_before,
    work_law.volume_after: _volume_after,
    work_law.pressure(work_law.volume): _pressure_expr
}).doit().simplify()

_work_from_law = law.rhs.subs({
    initial_volume: _volume_before,
    final_volume: _volume_after,
})

assert expr_equals(_work_expr, _work_from_law)


@validate_input(
    amount_of_substance_=amount_of_substance,
    start_volume_=initial_volume,
    final_volume_=final_volume,
    temperature_=temperature)
@validate_output(work)
def calculate_work(amount_of_substance_: Quantity, start_volume_: Quantity,
    final_volume_: Quantity, temperature_: Quantity) -> Quantity:
    solved = solve(law, work, dict=True)[0][work]
    result_expr = solved.subs({
        amount_of_substance: amount_of_substance_,
        initial_volume: start_volume_,
        final_volume: final_volume_,
        temperature: temperature_
    })
    return Quantity(result_expr)
