"""
Total energy transfer is zero in adiabatically isolated system
=============================================================

In an isolated system the total energy transferred between its parts is zero.
This is a direct consequence of the first law of thermodynamics.
"""

from sympy import (Eq, Idx, solve)
from symplyphysics import (units, Quantity, validate_input, validate_output, SymbolIndexed,
    SumIndexed, global_index)

amount_of_energy = SymbolIndexed("amount_of_energy", units.energy)
r"""
Amount of energy transferred between parts of the system.

Symbol:
    :code:`E_i`

Latex:
    :math:`E_i`
"""

law = Eq(SumIndexed(amount_of_energy[global_index], global_index), 0)
r"""
:code:`Sum(E_i, i) = 0`

Latex:
    .. math::
        \sum_i E_i = 0
"""


@validate_input(amounts_energy_=amount_of_energy)
@validate_output(units.energy)
def calculate_amount_energy(amounts_energy_: list[Quantity]) -> Quantity:
    local_index = Idx("index_local", (1, len(amounts_energy_) + 1))
    amounts_energy_law = law.subs(global_index, local_index)
    amounts_energy_law = amounts_energy_law.doit()
    unknown_amount_energy = amount_of_energy[len(amounts_energy_) + 1]
    solved = solve(amounts_energy_law, unknown_amount_energy, dict=True)[0][unknown_amount_energy]
    for i, v in enumerate(amounts_energy_):
        solved = solved.subs(amount_of_energy[i + 1], v)
    return Quantity(solved)
