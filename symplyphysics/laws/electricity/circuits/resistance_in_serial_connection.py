"""
Resistivity of serial resistors
===============================

The total resistance of the circuit whose components are connected in series is the sum
of the resistances of individual components.

**Conditions:**

#. Applies to direct current circuits.
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

total_resistance = symbols.resistance
"""
Total resistance of the circuit.
"""

resistance = SymbolIndexedNew("R_i", units.impedance)
r"""
Resistance of the :math:`i`-th component.
"""

law = Eq(total_resistance, SumIndexed(resistance[global_index], global_index))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(resistances_=resistance)
@validate_output(total_resistance)
def calculate_serial_resistance(resistances_: list[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(resistances_)))
    resistances_law = law.subs(global_index, local_index)
    resistances_law = resistances_law.doit()
    solved = solve(resistances_law, total_resistance, dict=True)[0][total_resistance]
    for i, v in enumerate(resistances_):
        solved = solved.subs(resistance[i + 1], v)
    return Quantity(solved)
