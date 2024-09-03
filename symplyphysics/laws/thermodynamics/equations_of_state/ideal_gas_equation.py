r"""
Ideal gas equation
==================

The ideal gas law, also known as the general gas equation, is an equation of state used to
describe a hypothetical ideal gas.

**Notation:**

#. :math:`R` is the molar gas constant.
"""

from sympy import (Eq, solve)
from symplyphysics import (symbols, units, Quantity, Symbol, validate_input, validate_output)

pressure = Symbol("pressure", units.pressure)
"""
Pressure inside the gas.

Symbol:
    :code:`p`
"""

volume = Symbol("volume", units.volume)
"""
Volume of the gas.

Symbol:
    :code:`V`
"""

amount_of_substance = Symbol("amount_of_substance", units.amount_of_substance)
"""
Amount of gas.

Symbol:
    :code:`n`
"""

temperature = symbols.thermodynamics.temperature
"""
:attr:`~symplyphysics.symbols.thermodynamics.temperature` of the gas.
"""

law = Eq(pressure * volume, amount_of_substance * units.molar_gas_constant * temperature)
r"""
:code:`p V = n * R * T`

Latex:
    .. math::
        p V = n R T
"""


@validate_input(volume_=volume, temperature_=temperature, mole_count_=amount_of_substance)
@validate_output(pressure)
def calculate_pressure(volume_: Quantity, temperature_: Quantity,
    mole_count_: Quantity) -> Quantity:
    solved = solve(law, pressure, dict=True)[0][pressure]
    result_expr = solved.subs({
        volume: volume_,
        temperature: temperature_,
        amount_of_substance: mole_count_
    })
    return Quantity(result_expr)
