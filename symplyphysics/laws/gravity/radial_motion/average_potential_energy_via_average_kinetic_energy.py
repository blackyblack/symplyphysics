r"""
Average potential energy via average kinetic energy
===================================================

The average potential energy of an orbiting planet is directly proportional to its average
kinetic energy. The average of both energies is taken w.r.t. time.

**Conditions:**

#. Works for elliptical orbits, i.e. the total energy :math:`E` of the planet is negative.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

average_potential_energy = Symbol("average_potential_energy", units.energy)
r"""
The potential energy of the planet averaged over time.

Symbol:
    :code:`avg(U)`

Latex:
    :math:`\langle U \rangle`
"""

average_kinetic_energy = Symbol("average_kinetic_energy", units.energy)
r"""
The kinetic energy of the planet averaged over time.

Symbol:
    :code:`avg(K)`

Latex:
    :math:`\langle K \rangle`
"""

law = Eq(average_potential_energy, -2 * average_kinetic_energy)
r"""
:code:`avg(U) = -2 * avg(K)`

Latex:
    .. math::
        \langle U \rangle = -2 \langle K \rangle
"""

# TODO Prove using the virial's theorem.


@validate_input(average_kinetic_energy_=average_kinetic_energy)
@validate_output(average_potential_energy)
def calculate_average_potential_energy(average_kinetic_energy_: Quantity) -> Quantity:
    result = law.rhs.subs(average_kinetic_energy, average_kinetic_energy_)
    return Quantity(result)
