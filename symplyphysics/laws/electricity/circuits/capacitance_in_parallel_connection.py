"""
Capacitance in parallel connection
==================================

The total capacitance of capacitors connected in parallel is the sum of the
capacitances of individual capacitors.
"""

from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, SymbolIndexed, SumIndexed, global_index)

total_capacitance = Symbol("total_capacitance", units.capacitance)
"""
Total capacitance.

Symbol:
    :code:`C`
"""

capacitance = SymbolIndexed("capacitance", units.capacitance)
r"""
Capacitance of :math:`i`-th capacitor.

Symbol:
    :code:`C_i`

Latex:
    :math:`C_i`
"""

law = Eq(total_capacitance, SumIndexed(capacitance[global_index], global_index))
r"""
:code:`C = Sum(C_i), i`

Latex:
    .. math::
        C = \sum_i C_i
"""


@validate_input(capacitances_=capacitance)
@validate_output(total_capacitance)
def calculate_parallel_capacitance(capacitances_: list[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(capacitances_)))
    capacitances_law = law.subs(global_index, local_index)
    capacitances_law = capacitances_law.doit()
    solved = solve(capacitances_law, total_capacitance, dict=True)[0][total_capacitance]
    for i, v in enumerate(capacitances_):
        solved = solved.subs(capacitance[i + 1], v)
    return Quantity(solved)
