"""
Current is voltage over resistance
==================================

Current flowing through a conductor is proportional to voltage applied to it and
inversely proportional to its resistance. This is also known as the **Ohm's law**.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

current = Symbol("current", units.current)
"""
Current flowing through the conductor.

Symbol:
    :code:`I`
"""

voltage = Symbol("voltage", units.voltage)
"""
Voltage applied to the conductor.

Symbol:
    :code:`V`
"""

resistance = Symbol("resistance", units.impedance)
r"""
Resistance of the conductor.

Symbol:
    :code:`R`
"""

law = Eq(current, voltage / resistance)
r"""
:code:`I = V / R`
"""


@validate_input(voltage_=voltage, resistance_=resistance)
@validate_output(current)
def calculate_current(voltage_: Quantity, resistance_: Quantity) -> Quantity:
    result_current_expr = solve(law, current, dict=True)[0][current]
    result_expr = result_current_expr.subs({voltage: voltage_, resistance: resistance_})
    return Quantity(result_expr)
