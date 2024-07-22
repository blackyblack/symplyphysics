"""
Total energy is negative average kinetic energy
===============================================

The total energy of an orbiting planet is equal to its negative kinetic energy averaged over time.

**Conditions:**

#. Works for elliptical (:math:`E < 0`) orbits.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

total_energy = Symbol("total_energy", units.energy)
"""
The total energy of the planet.

Symbol:
    E
"""

average_kinetic_energy = Symbol("average_kinetic_energy", units.energy)
"""
The kinetic energy of the planet averaged over time.

Symbol:
    K
"""

law = Eq(total_energy, -1 * average_kinetic_energy)
"""
E = -1 * K

Latex:
    .. math::
        E = -K
"""


@validate_input(average_kinetic_energy_=average_kinetic_energy)
@validate_output(total_energy)
def calculate_total_energy(average_kinetic_energy_: Quantity) -> Quantity:
    result = law.rhs.subs(average_kinetic_energy, average_kinetic_energy_)
    return Quantity(result)
