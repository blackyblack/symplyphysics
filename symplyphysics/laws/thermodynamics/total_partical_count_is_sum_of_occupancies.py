r"""
Total particle count is sum of occupancies
==========================================

Maxwell—Boltzmann, Fermi—Dirac and Bose—Einstein distributions describe the average number of particles
in some energy level or particle state. They are normalized by the condition that the sum of the number
of particles in each energy level or particle state should be equal to the total number of particles
in the system.

**Notes:**

#. This law can be used to express chemical potential :math:`\mu` as a function of temperature :math:`T`
   and total particle count :math:`N` for Fermi—Dirac and Bose—Einstein distributions.

**Links:**

#. `Wikipedia, see normalization condition (second formula) <https://en.wikipedia.org/wiki/Fermi%E2%80%93Dirac_statistics#Fermi%E2%80%93Dirac_distribution>`__.
"""

from typing import Sequence
from sympy import Eq, Idx
from symplyphysics import (
    dimensionless,
    validate_input,
    validate_output,
    global_index,
    SymbolIndexedNew,
    SumIndexed,
    symbols,
)

total_particle_count = symbols.particle_count
"""
Total :symbols:`particle_count` of the system.
"""

occupancy = SymbolIndexedNew("N", dimension=dimensionless)
"""
Occupancy of energy level or particle state of index :math:`i`.
"""

law = Eq(total_particle_count, SumIndexed(occupancy[global_index], global_index))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(occupancies_=occupancy)
@validate_output(total_particle_count)
def calculate_total_particle_count(occupancies_: Sequence[float]) -> int:
    local_index_ = Idx("local_index_", (1, len(occupancies_)))
    result = law.rhs.subs(global_index, local_index_).doit()
    for idx_, count_ in enumerate(occupancies_, 1):
        result = result.subs(occupancy[idx_], count_)
    return int(result)
