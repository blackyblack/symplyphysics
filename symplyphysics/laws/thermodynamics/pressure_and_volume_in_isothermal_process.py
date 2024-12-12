"""
Pressure and volume in isothermal process
=========================================

A thermodynamic process is isothermal when the temperature of the system stays the same.
**Boyle's law** states that the product of pressure and volume is constant during an isothermal
process.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Ideal_gas_law#math_1>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.laws.thermodynamics.equations_of_state import (
    ideal_gas_equation as thermodynamics_law,
)

initial_pressure = clone_as_symbol(symbols.pressure, subscript="0")
"""
Initial :symbols:`pressure` inside the system.
"""

final_pressure = clone_as_symbol(symbols.pressure, subscript="1")
"""
Final :symbols:`pressure` inside the system.
"""

initial_volume = clone_as_symbol(symbols.volume, subscript="0")
"""
Initial :symbols:`volume` of the system.
"""

final_volume = clone_as_symbol(symbols.volume, subscript="1")
"""
Final :symbols:`volume` of the system.
"""

law = Eq(initial_pressure * initial_volume, final_pressure * final_volume)
"""
:laws:symbol::

:laws:latex::
"""

## Derive the same law from the general ideal gas law

_temperature_start = clone_as_symbol(symbols.temperature, subscript="0")
_temperature_end = clone_as_symbol(symbols.temperature, subscript="1")

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
