"""
Power via voltage and current
=============================

Electric power can be expressed using current flowing through a
conductor and voltage applied.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols

power = symbols.power
"""
Electric :symbols:`power`.
"""

current = symbols.current
"""
Electric ;symbols:`current`.
"""

voltage = symbols.voltage
"""
:symbols:`voltage` applied.
"""

law = Eq(power, current * voltage)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(current_=current, voltage_=voltage)
@validate_output(power)
def calculate_power(current_: Quantity, voltage_: Quantity) -> Quantity:
    result_power_expr = solve(law, power, dict=True)[0][power]
    result_expr = result_power_expr.subs({current: current_, voltage: voltage_})
    return Quantity(result_expr)
