"""
Rotational inertia is additive
==============================

For a system composed of several parts, its total rotational inertia is the sum of the rotational
inertia of each part of the system.

**Conditions:**

#. The rotational inertia is calculated for the same axis for all parts of the system.
"""

from typing import Sequence
from sympy import Eq, Idx, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    SymbolIndexed,
    SumIndexed,
    global_index,
)

total_rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
"""
Total rotational inertia of the system.

Symbol:
    :code:`I`
"""

rotational_inertia = SymbolIndexed("rotational_inertia", units.mass * units.length**2)
r"""
Rotational inertia of the :math:`k`-th part of the system.

Symbol:
    :code:`I_k`

Latex:
    :math:`I_k`
"""

law = Eq(total_rotational_inertia, SumIndexed(rotational_inertia[global_index], global_index))
r"""
:code:`I = Sum(I_k, k)`

Latex:
    .. math::
        I = \sum_k I_k
"""


@validate_input(rotational_inertias_=rotational_inertia)
@validate_output(total_rotational_inertia)
def calculate_rotational_inertia(rotational_inertias_: Sequence[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(rotational_inertias_)))
    rotational_inertia_law = law.subs(global_index, local_index)
    rotational_inertia_law = rotational_inertia_law.doit()
    solved = solve(rotational_inertia_law, total_rotational_inertia,
        dict=True)[0][total_rotational_inertia]
    for i, v in enumerate(rotational_inertias_):
        solved = solved.subs(rotational_inertia[i + 1], v)
    return Quantity(solved)
