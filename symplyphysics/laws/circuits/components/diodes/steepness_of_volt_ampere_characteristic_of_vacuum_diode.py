"""
Steepness of volt-ampere characteristic of vacuum diode
=======================================================

The steepness of the volt-ampere characteristic is generally equal to the voltage
derivative of the current. For a vacuum diode, it can be calculated if the diode
parameter and the voltage between the anode and cathode are known.

..
    TODO: find link
"""

from sympy import Eq, Rational, solve, sqrt
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
)

steepness = Symbol("S", units.current / units.voltage)
"""
Steepness of the volt-ampere characteristic of the vacuum diode.
"""

diode_constant = symbols.diode_constant
"""
:symbols:`diode_constant`.
"""

anode_voltage = symbols.voltage
"""
:symbols:`voltage` between cathode and anode.
"""

law = Eq(steepness, Rational(3, 2) * diode_constant * sqrt(anode_voltage))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(diode_constant_=diode_constant, voltage_=anode_voltage)
@validate_output(steepness)
def calculate_steepness(diode_constant_: Quantity, voltage_: Quantity) -> Quantity:
    result_expr = solve(law, steepness, dict=True)[0][steepness]
    result_expr = result_expr.subs({
        diode_constant: diode_constant_,
        anode_voltage: voltage_,
    })
    return Quantity(result_expr)
