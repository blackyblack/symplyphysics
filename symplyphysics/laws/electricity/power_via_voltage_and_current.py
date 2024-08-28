"""
Power via voltage and current
=============================

Electric power can be expressed using current flowing through a
conductor and voltage applied.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

power = Symbol("power", units.power)
"""
Electric power.

Symbol:
    :code:`P`
"""

current = Symbol("current", units.current)
"""
Electric current.

Symbol:
    :code:`I`
"""

voltage = Symbol("voltage", units.voltage)
"""
Voltage applied.

Symbol:
    :code:`V`
"""

law = Eq(power, current * voltage)
r"""
:code:`P = I * V`

Latex:
    .. math::
        P = I V
"""


@validate_input(current_=current, voltage_=voltage)
@validate_output(power)
def calculate_power(current_: Quantity, voltage_: Quantity) -> Quantity:
    result_power_expr = solve(law, power, dict=True)[0][power]
    result_expr = result_power_expr.subs({current: current_, voltage: voltage_})
    return Quantity(result_expr)
