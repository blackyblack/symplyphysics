"""
Admittance in parallel connection
=================================

The total admittance of the circuit whose components are connected in parallel is the sum
of the admittances of individual components.
"""

from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, global_index, SumIndexed, SymbolIndexed)

total_admittance = Symbol("total_admittance", units.conductance)
"""
Total admittance of the circuit.

Symbol:
    :code:`Y`
"""

admittance = SymbolIndexed("admittance", units.conductance)
r"""
Admittance of :math:`i`-th circuit.

Symbol:
    :code:`Y_i`

Latex:
    :math:`Y_i`
"""

law = Eq(total_admittance, SumIndexed(admittance[global_index], global_index))
r"""
:code:`Y = Sum(Y_i, i)`

Latex:
    .. math::
        Y = \sum_i Y_i
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
