"""
Power via voltage and resistance
================================

Power can be found using voltage and resistance.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

power = Symbol("power", units.power)
"""
Power.

Symbol:
    :code:`P`
"""

voltage = Symbol("voltage", units.voltage)
"""
Voltage.

Symbol:
    :code:`V`
"""

resistance = Symbol("resistance", units.impedance)
"""
Resistance.

Symbol:
    :code:`R`
"""

law = Eq(power, voltage**2 / resistance)
r"""
:code:`P = V^2 / R`

Latex:
    .. math::
        P = \frac{V^2}{R}
"""


@validate_input(
    voltage_=voltage,
    resistance_=resistance,
)
@validate_output(power)
def calculate_power(
    voltage_: Quantity,
    resistance_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        voltage: voltage_,
        resistance: resistance_,
    })
    return Quantity(result)
