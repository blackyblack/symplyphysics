"""
Gas mixture pressure from partial pressures
===========================================

The pressure of a mixture of gases is equal to the sum of their partial pressures. *Partial pressure*
of a gas is the pressure that the gas would exert on the walls of a vessel when it is the only gas present
in the vessel.
"""

from typing import Sequence
from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, SymbolIndexed, SumIndexed, global_index)

total_pressure = Symbol("total_pressure", units.pressure)
"""
Total pressure inside the system.

Symbol:
    :code:`p`
"""

partial_pressure = SymbolIndexed("partial_pressure", units.pressure)
r"""
Partial pressure of the :math:`i`-th gas component.

Symbol:
    :code:`p_i`

Latex:
    :math:`p_i`
"""

law = Eq(total_pressure, SumIndexed(partial_pressure[global_index], global_index))
r"""
:code:`p = Sum(p_i, i)`

Latex:
    .. math::
        p = \sum_i p_i
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
