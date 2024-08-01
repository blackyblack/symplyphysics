"""
Chemical potential is Gibbs energy per particle
===============================================

Chemical potential of a thermodynamic system can be expressed as the Gibbs energy of the system
per particle. Therefore, chemical potential is an intensive physical quantity, whereas Gibbs energy
and particle count are extensive.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

chemical_potential = Symbol("chemical_potential", units.energy)
r"""
Chemical potential of the system.

Symbol:
    :code:`mu`

Latex:
    :math:`\mu`
"""

gibbs_energy = Symbol("gibbs_energy", units.energy)
"""
Gibbs energy of the system.

Symbol:
    :code:`G`
"""

particle_count = Symbol("particle_count")
"""
Number of particles in the system.

Symbol:
    :code:`N`
"""

law = Eq(chemical_potential, gibbs_energy / particle_count)
r"""
:code:`mu = G / N`

Latex:
    .. math::
        \mu = \frac{G}{N}
"""

# TODO: derive law from the extensive property of chemical potential


@validate_input(
    gibbs_energy_=gibbs_energy,
    particle_count_=particle_count,
)
@validate_output(chemical_potential)
def calculate_chemical_potential(
    gibbs_energy_: Quantity,
    particle_count_: int,
) -> Quantity:
    result = law.rhs.subs({
        gibbs_energy: gibbs_energy_,
        particle_count: particle_count_,
    })
    return Quantity(result)
