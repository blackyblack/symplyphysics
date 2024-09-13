"""
Sum of currents through junction is zero
========================================

Electrical charge is neither created nor accumulated in electrical circuits.
Alternatively, one can think that the total current flowing into an electrical
junction must be equal to the total current flowing out of it. This law is also
known as the **first Kirchhoff's law**, or the **Kirchhoff's junction rule**.

**Notes:**

#. Current flowing into the junction and current flowing out of the junction
   are of opposite signs.
"""

from typing import Sequence
from sympy import Eq, Idx, solve
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    SymbolIndexedNew,
    SumIndexed,
    global_index,
)

current = SymbolIndexedNew("I_k", units.current)
r"""
:math:`k`-th current flowing through the node.
"""

law = Eq(SumIndexed(current[global_index], global_index), 0)
"""
:laws:symbol::

:laws:latex::
"""

# TODO Derive from law of conservation of charge


@validate_input(currents_=current)
@validate_output(units.current)
def calculate_current_from_array(currents_: Sequence[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(currents_) + 1))
    currents_law = law.subs(global_index, local_index)
    currents_law = currents_law.doit()
    unknown_current = current[len(currents_) + 1]
    solved = solve(currents_law, unknown_current, dict=True)[0][unknown_current]
    for i, v in enumerate(currents_):
        solved = solved.subs(current[i + 1], v)
    return Quantity(solved)
