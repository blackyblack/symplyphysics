"""
Volume and temperature in isobaric process
==========================================

A thermodynamic process is called **isobaric** when the pressure inside the system stays
constant. Also called **Charles's law**, it states that the volume of the gas scales
by the same amount as the temperature during an isobaric process.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Ideal_gas_law#math_2>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
)
from symplyphysics.laws.thermodynamics.equations_of_state import (
    ideal_gas_equation as thermodynamics_law,)

initial_temperature = clone_as_symbol(symbols.temperature, subscript="0")
"""
Initial :symbols:`temperature` of the system.
"""

final_temperature = clone_as_symbol(symbols.temperature, subscript="1")
"""
Final :symbols:`temperature` of the system.
"""

initial_volume = clone_as_symbol(symbols.volume, subscript="0")
"""
Initial :symbols:`volume` of the system.
"""

final_volume = clone_as_symbol(symbols.volume, subscript="1")
"""
Final :symbols:`volume` of the system.
"""

law = Eq(initial_volume / final_volume, initial_temperature / final_temperature)
"""
:laws:symbol::

:laws:latex::
"""

## Derive the same law from the general ideal gas law

_pressure_start = clone_as_symbol(symbols.pressure, subscript="0")
_pressure_end = clone_as_symbol(symbols.pressure, subscript="1")

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
