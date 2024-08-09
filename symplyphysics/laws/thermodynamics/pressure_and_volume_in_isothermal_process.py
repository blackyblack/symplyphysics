"""
Pressure and volume in isothermal process
=========================================

A thermodynamic process is isothermal when the temperature of the system stays the same.
Boyle's law states that the product of pressure and volume is constant during an isothermal
process.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation as thermodynamics_law

initial_pressure = Symbol("initial_pressure", units.pressure)
r"""
Initial pressure inside the system.

Symbol:
    :code:`p_0`

Latex:
    :math:`p_0`
"""

final_pressure = Symbol("final_pressure", units.pressure)
r"""
Final pressure inside the system.

Symbol:
    :code:`p_1`

Latex:
    :math:`p_1`
"""

initial_volume = Symbol("initial_volume", units.volume)
r"""
Initial volume of the system.

Symbol:
    :code:`V_0`

Latex:
    :math:`V_0`
"""

final_volume = Symbol("final_volume", units.volume)
r"""
Final volume of the system.

Symbol:
    :code:`V_1`

Latex:
    :math:`V_1`
"""

law = Eq(initial_pressure * initial_volume, final_pressure * final_volume)
r"""
:code:`p_0 V_0 = p_1 V_1`

Latex:
    .. math::
        p_0 V_0 = p_1 V_1
"""

## Derive the same law from the general ideal gas law

_temperature_start = Symbol("temperature_start", units.temperature)
_temperature_end = Symbol("temperature_end", units.temperature)

_isothermal_condition = Eq(_temperature_start, _temperature_end)

_eq_start = thermodynamics_law.law.subs({
    thermodynamics_law.temperature: _temperature_start,
    thermodynamics_law.volume: initial_volume,
    thermodynamics_law.pressure: initial_pressure
})

_eq_end = thermodynamics_law.law.subs({
    thermodynamics_law.temperature: _temperature_end,
    thermodynamics_law.volume: final_volume,
    thermodynamics_law.pressure: final_pressure
})

_derived_law = [_eq_start, _eq_end, _isothermal_condition]

## Check the equivalence of 'law' and '_derived_law'
_derived_pressure_end = solve(_derived_law, (_temperature_start, _temperature_end, final_pressure),
    dict=True)[0][final_pressure]
assert solve(law, final_pressure, dict=True)[0][final_pressure] == _derived_pressure_end


@validate_input(pressure_start_=initial_pressure,
    pressure_end_=final_pressure,
    volume_start_=initial_volume)
@validate_output(final_volume)
def calculate_volume(pressure_start_: Quantity, volume_start_: Quantity,
    pressure_end_: Quantity) -> Quantity:
    solved = solve(law, final_volume, dict=True)[0][final_volume]
    result_expr = solved.subs({
        initial_pressure: pressure_start_,
        initial_volume: volume_start_,
        final_pressure: pressure_end_
    })
    return Quantity(result_expr)
