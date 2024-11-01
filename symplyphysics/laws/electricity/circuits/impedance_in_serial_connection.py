"""
Impedance in serial connection
==============================

The total impedance of a circuit whose components are connected in series is the sum
of the impedances of individual components.

**Conditions:**

#. Components are connected in series.
"""

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

total_impedance = symbols.electrical_impedance
"""
Total :symbols:`electrical_impedance` of the circuit.
"""

impedance = clone_as_indexed(symbols.electrical_impedance,
    display_symbol="Z[i]",
    display_latex="Z_i")
r"""
:symbols:`electrical_impedance` of the :math:`i`-th component.
"""

law = Eq(total_impedance, SumIndexed(impedance[global_index], global_index))
"""
:laws:symbol::

:laws:latex::
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
