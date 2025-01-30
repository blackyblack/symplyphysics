"""
Gas mixture pressure from partial pressures
===========================================

The pressure of a mixture of gases is equal to the sum of their partial pressures. *Partial pressure*
of a gas is the pressure that the gas would exert on the walls of a vessel when it is the only gas present
in the vessel.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Dalton%27s_law>`__.
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

total_pressure = symbols.pressure
"""
Total :symbols:`pressure` inside the system.

Symbol:
    :code:`p`
"""

partial_pressure = clone_as_indexed(symbols.pressure)
"""
Partial :symbols:`pressure` of the :math:`i`-th gas component.
"""

law = Eq(total_pressure, SumIndexed(partial_pressure[global_index], global_index))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(partial_pressures_=partial_pressure)
@validate_output(total_pressure)
def calculate_total_pressure(partial_pressures_: Sequence[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(partial_pressures_)))
    partial_pressures_law = law.subs(global_index, local_index)
    partial_pressures_law = partial_pressures_law.doit()
    solved = solve(partial_pressures_law, total_pressure, dict=True)[0][total_pressure]
    for i, v in enumerate(partial_pressures_):
        solved = solved.subs(partial_pressure[i + 1], v)
    return Quantity(solved)
