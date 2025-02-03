"""
Internal resistance of vacuum diode
===================================

The internal resistance is generally equal to the current derivative of the voltage. For
a vacuum diode, it can be calculated if the diode parameter and the voltage between the
anode and cathode are known.

..
    TODO: find link
"""

from sympy import Eq, solve, sqrt
from symplyphysics import Quantity, validate_input, validate_output, symbols

internal_resistance = symbols.electrical_resistance
"""
Internal :symbols:`electrical_resistance` of the vacuum diode.
"""

diode_constant = symbols.diode_constant
"""
:symbols:`diode_constant`.
"""

voltage = symbols.voltage
"""
:symbols:`voltage` between cathode and anode.
"""

law = Eq(internal_resistance, 2 / (3 * diode_constant * sqrt(voltage)))


@validate_input(diode_constant_=diode_constant, voltage_=voltage)
@validate_output(internal_resistance)
def calculate_internal_resistance(diode_constant_: Quantity, voltage_: Quantity) -> Quantity:
    result_expr = solve(law, internal_resistance, dict=True)[0][internal_resistance]
    result_expr = result_expr.subs({
        diode_constant: diode_constant_,
        voltage: voltage_,
    })
    return Quantity(result_expr)
