"""
Admittance in parallel connection
=================================

The total admittance of the circuit whose components are connected in parallel is the sum
of the admittances of individual components.

**Conditions:**

#. Components are connected in parallel.
"""

from sympy import Eq, Idx, solve
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    global_index,
    SumIndexed,
    symbols,
)
from symplyphysics.core.symbols.symbols import clone_as_indexed

total_admittance = symbols.admittance
"""
Total :symbols:`admittance` of the circuit.
"""

admittance = clone_as_indexed(symbols.admittance)
"""
:symbols:`admittance` of :math:`i`-th circuit.
"""

law = Eq(total_admittance, SumIndexed(admittance[global_index], global_index))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(admittances_=admittance)
@validate_output(units.conductance)
def calculate_parallel_admittance(admittances_: list[Quantity]) -> Quantity:
    local_index = Idx("local_index", (1, len(admittances_)))
    admittances_law = law.subs(global_index, local_index)
    admittances_law = admittances_law.doit()
    solved = solve(admittances_law, total_admittance, dict=True)[0][total_admittance]
    for i, v in enumerate(admittances_):
        solved = solved.subs(admittance[i + 1], v)
    return Quantity(solved)
