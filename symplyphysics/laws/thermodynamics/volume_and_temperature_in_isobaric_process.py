"""
Volume and temperature in isobaric process
==========================================

A thermodynamic process is called isobaric when the pressure inside the system stays
constant. Also called Guy-Lussac's law, it states that the volume of the gas scales
by the same amount as the temperature during an isobaric process.
"""

from sympy import (Eq, solve)
from symplyphysics import (clone_as_symbol, symbols, units, Quantity, Symbol, validate_input,
    validate_output)
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation as thermodynamics_law

initial_temperature = clone_as_symbol(symbols.temperature,
    display_symbol="T0",
    display_latex="T_0")
"""
Initial :attr:`~symplyphysics.symbols.temperature` of the system.
"""

final_temperature = clone_as_symbol(symbols.temperature,
    display_symbol="T1",
    display_latex="T_1")
"""
Final :attr:`~symplyphysics.symbols.temperature` of the system.
"""

initial_volume = Symbol("initial_volume", units.volume)
"""
Initial volume of the system.

Symbol: 
    :code:`V0`

Latex:
    :math:`V_0`
"""

final_volume = Symbol("final_volume", units.volume)
"""
Final volume of the system.

Symbol: 
    :code:`V1`

Latex:
    :math:`V_1`
"""

law = Eq(initial_volume / final_volume, initial_temperature / final_temperature)
r"""
:code:`V0 / V1 = T0 / T1`

Latex:
    .. math::
        \frac{V_0}{V_1} = \frac{T_0}{T_1}
"""

## Derive the same law from the general ideal gas law

_pressure_start = Symbol("_pressure_start", units.pressure)
_pressure_end = Symbol("_pressure_end", units.pressure)

_isobaric_condition = Eq(_pressure_start, _pressure_end)

_eq_start = thermodynamics_law.law.subs({
    thermodynamics_law.temperature: initial_temperature,
    thermodynamics_law.volume: initial_volume,
    thermodynamics_law.pressure: _pressure_start
})

_eq_end = thermodynamics_law.law.subs({
    thermodynamics_law.temperature: final_temperature,
    thermodynamics_law.volume: final_volume,
    thermodynamics_law.pressure: _pressure_end
})

_derived_law = [_eq_start, _eq_end, _isobaric_condition]

## Check the equivalence of 'law' and '_derived_law'
_derived_temperature_end = solve(_derived_law, (_pressure_start, _pressure_end, final_temperature),
    dict=True)[0][final_temperature]
assert solve(law, final_temperature, dict=True)[0][final_temperature] == _derived_temperature_end


@validate_input(temperature_start_=initial_temperature,
    temperature_end_=final_temperature,
    volume_start_=initial_volume)
@validate_output(final_volume)
def calculate_volume(temperature_start_: Quantity, volume_start_: Quantity,
    temperature_end_: Quantity) -> Quantity:
    solved = solve(law, final_volume, dict=True)[0][final_volume]
    result_expr = solved.subs({
        initial_temperature: temperature_start_,
        initial_volume: volume_start_,
        final_temperature: temperature_end_
    })
    return Quantity(result_expr)
