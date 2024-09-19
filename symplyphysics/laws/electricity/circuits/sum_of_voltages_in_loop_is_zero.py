"""
Sum of voltages in loop is zero
===============================

The directed sum of the potential differences, or voltages, around any closed loop
is zero. Directed sum implies that the sign of the voltages must be taken into account.
This law is also known as **Kirchhoff's second law**, or **Kirchhoff's loop rule**.
"""

from typing import Sequence
from sympy import Eq, Idx, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    SumIndexed,
    global_index,
    symbols,
)
from symplyphysics.core.symbols.symbols import clone_as_indexed

voltage = clone_as_indexed(symbols.voltage, display_symbol="V[i]", display_latex="V_i")
r"""
:math:`i`-th voltage.
"""

law = Eq(SumIndexed(voltage[global_index], global_index), 0)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(voltages_=voltage)
@validate_output(voltage)
def calculate_voltage(voltages_: Sequence[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(voltages_) + 1))
    voltages_law = law.subs(global_index, local_index)
    voltages_law = voltages_law.doit()
    unknown_voltage = voltage[len(voltages_) + 1]
    solved = solve(voltages_law, unknown_voltage, dict=True)[0][unknown_voltage]
    for i, v in enumerate(voltages_):
        solved = solved.subs(voltage[i + 1], v)
    return Quantity(solved)
