"""
Rotational inertia is additive
==============================

For a system composed of several parts, its total rotational inertia is the sum of the rotational
inertia of each part of the system.

**Conditions:**

#. The rotational inertia is calculated for the same axis for all parts of the system.

**Links:**

#. `Wikipedia, see second paragraph <https://en.wikipedia.org/wiki/Moment_of_inertia>`__.
"""

from typing import Sequence
from sympy import Eq, Idx, solve
from symplyphysics import Quantity, validate_input, validate_output, SumIndexed, symbols
from symplyphysics.core.symbols.symbols import clone_as_indexed

total_rotational_inertia = symbols.rotational_inertia
"""
Total :symbols:`rotational_inertia` of the system.
"""

_index = Idx("k")

rotational_inertia = clone_as_indexed(symbols.rotational_inertia, _index)
"""
:symbols:`rotational_inertia` of the :math:`k`-th part of the system.
"""

law = Eq(total_rotational_inertia, SumIndexed(rotational_inertia[_index], _index))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(rotational_inertias_=rotational_inertia)
@validate_output(total_rotational_inertia)
def calculate_rotational_inertia(rotational_inertias_: Sequence[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(rotational_inertias_)))
    rotational_inertia_law = law.subs(_index, local_index)
    rotational_inertia_law = rotational_inertia_law.doit()
    solved = solve(rotational_inertia_law, total_rotational_inertia,
        dict=True)[0][total_rotational_inertia]
    for i, v in enumerate(rotational_inertias_):
        solved = solved.subs(rotational_inertia[i + 1], v)
    return Quantity(solved)
