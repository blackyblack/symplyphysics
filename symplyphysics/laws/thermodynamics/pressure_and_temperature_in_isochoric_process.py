"""
Pressure and temperature in isochoric process
=============================================

A thermodynamic process is isochoric when the volume of the body stays constant during it.
During an isochoric process, the pressure of the gas scales by the same amount as the temperature
of the gas.

**Condition:**

#. Applies to ideal gases.
"""

from sympy import (Eq, solve)
from symplyphysics import (clone_symbol, symbols, units, Quantity, Symbol,
    validate_input, validate_output)
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation as thermodynamics_law

initial_pressure = Symbol("initial_pressure", units.pressure)
r"""
Initial pressure of the gas.

Symbol:
    :code:`p_0`

Latex:
    :math:`p_0`
"""

final_pressure = Symbol("final_pressure", units.pressure)
r"""
Final pressure of the gas.

Symbol:
    :code:`p_1`

Latex:
    :math:`p_1`
"""

initial_temperature = clone_symbol(symbols.thermodynamics.temperature, "initial_temperature")
r"""
Initial :attr:`~symplyphysics.symbols.thermodynamics.temperature` of the gas.

Symbol:
    :code:`T_0`

Latex:
    :math:`T_0`
"""

final_temperature = clone_symbol(symbols.thermodynamics.temperature, "final_temperature")
r"""
Final :attr:`~symplyphysics.symbols.thermodynamics.temperature` of the gas.

Symbol:
    :code:`T_1`

Latex:
    :math:`T_1`
"""

law = Eq(initial_pressure / final_pressure, initial_temperature / final_temperature)
r"""
:code:`p_0 / p_1 = T_0 / T_1`

Latex:
    .. math::
        \frac{p_0}{p_1} = \frac{T_0}{T_1}
"""

## Derive the same law from the general ideal gas law

_volume_start = Symbol("volume_start", units.volume)
_volume_end = Symbol("volume_end", units.volume)

_isochoric_condition = Eq(_volume_start, _volume_end)

_eq_start = thermodynamics_law.law.subs({
    thermodynamics_law.temperature: initial_temperature,
    thermodynamics_law.volume: _volume_start,
    thermodynamics_law.pressure: initial_pressure
})

_eq_end = thermodynamics_law.law.subs({
    thermodynamics_law.temperature: final_temperature,
    thermodynamics_law.volume: _volume_end,
    thermodynamics_law.pressure: final_pressure
})

_derived_law = [_eq_start, _eq_end, _isochoric_condition]

## Check the equivalence of 'law' and '_derived_law'
_derived_pressure_end = solve(_derived_law, (_volume_start, _volume_end, final_pressure),
    dict=True)[0][final_pressure]

assert solve(law, final_pressure, dict=True)[0][final_pressure] == _derived_pressure_end


@validate_input(temperature_start_=initial_temperature,
    pressure_start_=initial_pressure,
    temperature_end_=final_temperature)
@validate_output(final_pressure)
def calculate_pressure(temperature_start_: Quantity, pressure_start_: Quantity,
    temperature_end_: Quantity) -> Quantity:
    solved = solve(law, final_pressure, dict=True)[0][final_pressure]
    result_expr = solved.subs({
        initial_pressure: pressure_start_,
        initial_temperature: temperature_start_,
        final_temperature: temperature_end_
    })
    return Quantity(result_expr)