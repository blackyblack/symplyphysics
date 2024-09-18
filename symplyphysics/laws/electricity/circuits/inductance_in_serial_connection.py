"""
Inductance in serial connection
===============================

The total inductance of the circuit whose components are connected in series is the sum
of the inductances of individual components.

**Conditions:**

#. Components are connected in series.
#. Inductors are not magnetically coupled.
"""

from sympy import Eq, Idx, solve
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    SymbolIndexedNew,
    SumIndexed,
    global_index,
    symbols,
)

total_inductance = symbols.inductance
"""
Total :symbols:`inductance` of the circuit.
"""

inductance = SymbolIndexedNew("L[i]", units.inductance, display_latex="L_i")
r"""
:symbols:`inductance` of the :math:`i`-th component.
"""

law = Eq(total_inductance, SumIndexed(inductance[global_index], global_index))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(inductances_=inductance)
@validate_output(total_inductance)
def calculate_serial_inductance(inductances_: list[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(inductances_)))
    inductances_law = law.subs(global_index, local_index)
    inductances_law = inductances_law.doit()
    solved = solve(inductances_law, total_inductance, dict=True)[0][total_inductance]
    for i, v in enumerate(inductances_):
        solved = solved.subs(inductance[i + 1], v)
    return Quantity(solved)
