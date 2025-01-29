"""
Total energy transfer is zero in adiabatically isolated system
=============================================================

In an isolated system the total energy transferred between its parts is zero.
This is a direct consequence of the first law of thermodynamics.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Isolated_system>`__.
"""

from sympy import Eq, Idx, solve
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    SumIndexed,
    global_index,
    symbols,
)
from symplyphysics.core.symbols.symbols import clone_as_indexed

energy = clone_as_indexed(symbols.energy)
"""
Amount of :symbols:`energy` transferred between parts of the system.
"""

law = Eq(SumIndexed(energy[global_index], global_index), 0)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(amounts_energy_=energy)
@validate_output(units.energy)
def calculate_amount_energy(amounts_energy_: list[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(amounts_energy_) + 1))
    amounts_energy_law = law.subs(global_index, local_index)
    amounts_energy_law = amounts_energy_law.doit()
    unknown_amount_energy = energy[len(amounts_energy_) + 1]
    solved = solve(amounts_energy_law, unknown_amount_energy, dict=True)[0][unknown_amount_energy]
    for i, v in enumerate(amounts_energy_):
        solved = solved.subs(energy[i + 1], v)
    return Quantity(solved)
