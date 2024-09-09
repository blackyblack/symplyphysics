"""
Impedance in serial connection
==============================

The total impedance of a circuit whose components are connected in series is the sum
of the impedances of individual components.

**Conditions:**

#. Components are connected in series.
"""

from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, SymbolIndexed, SumIndexed, global_index)

total_impedance = Symbol("total_impedance", units.impedance)
"""
Total impedance of the circuit.

Symbol:
    :code:`Z`
"""

impedance = SymbolIndexed("impedance", units.impedance)
r"""
Impedance of the :math:`i`-th component.

Symbol:
    :code:`Z_i`

Latex:
    :math:`Z_i`
"""

law = Eq(total_impedance, SumIndexed(impedance[global_index], global_index))
r"""
:code:`Z = Sum(Z_i, i)`

Latex:
    .. math::
        Z = \sum_i Z_i
"""


@validate_input(impedances_=impedance)
@validate_output(total_impedance)
def calculate_serial_impedance(impedances_: list[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(impedances_)))
    impedances_law = law.subs(global_index, local_index)
    impedances_law = impedances_law.doit()
    solved = solve(impedances_law, total_impedance, dict=True)[0][total_impedance]
    for i, v in enumerate(impedances_):
        solved = solved.subs(impedance[i + 1], v)
    return Quantity(solved)
