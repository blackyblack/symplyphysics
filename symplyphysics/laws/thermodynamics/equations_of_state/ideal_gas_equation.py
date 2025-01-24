r"""
Ideal gas equation
==================

The ideal gas law, also known as the general gas equation, is an equation of state used to
describe a hypothetical ideal gas.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Ideal_gas_law>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
    quantities,
)

pressure = symbols.pressure
"""
:symbols:`pressure` inside the gas.
"""

volume = symbols.volume
"""
:symbols:`volume` of the gas.
"""

amount_of_substance = symbols.amount_of_substance
"""
:symbols:`amount_of_substance` of the gas.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the gas.
"""

law = Eq(pressure * volume, amount_of_substance * quantities.molar_gas_constant * temperature)
"""
:laws:symbol::

:laws:latex::
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
