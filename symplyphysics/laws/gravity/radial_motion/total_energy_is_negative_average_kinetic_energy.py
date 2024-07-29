"""
Total energy is negative average kinetic energy
===============================================

The total energy of an orbiting planet is equal to its negative kinetic energy averaged over time.

**Conditions:**

#. Works for elliptical (:math:`E < 0`) orbits.

**Notes:**

#. :code:`<A>` denotes the average of the quantity :code:`A`.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

total_mechanical_energy = Symbol("total_mechanical_energy", units.energy)
"""
The total mechanical energy of the planet.

Symbol:
    :code:`E`
"""

average_kinetic_energy = Symbol("average_kinetic_energy", units.energy)
r"""
The kinetic energy of the planet averaged over time.

Symbol:
    :code:`<K>`

Latex:
    :math:`\langle K \rangle`
"""

law = Eq(total_mechanical_energy, -1 * average_kinetic_energy)
r"""
:code:`E = -1 * <K>`

Latex:
    .. math::
        E = - \langle K \rangle
"""


@validate_input(average_kinetic_energy_=average_kinetic_energy)
@validate_output(total_mechanical_energy)
def calculate_total_energy(average_kinetic_energy_: Quantity) -> Quantity:
    result = law.rhs.subs(average_kinetic_energy, average_kinetic_energy_)
    return Quantity(result)
