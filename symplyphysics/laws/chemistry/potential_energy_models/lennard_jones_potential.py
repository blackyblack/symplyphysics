"""
Lennard-Jones potential
=======================

The Lennard-Jones (LJ) potential is an intermolecular pair potential. It models soft repulsive and
attractive (van der Waals) interactions, and it describes electronically neutral atoms or molecules.
It is a simple yet realistic model to describe intermolecular interactions: two particles repel each
other at a very close distance, attract each other at moderate distance, and do not interact at infinite
distance.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

potential = Symbol("potential", units.energy)
"""
Potential energy of the configuration.

Symbol:
    :code:`U`
"""

dispersion_energy = Symbol("dispersion_energy", units.energy)
r"""
Depth of the potential well, also referred to as "dispersion energy".

Symbol:
    :code:`e`

Latex:
    :math:`\varepsilon`
"""

particle_size = Symbol("particle_size", units.length)
r"""
Distance at which potential is zero, also referred to as "particle size".

Symbol:
    :code:`sigma`

Latex:
    :math:`\sigma`
"""

distance = Symbol("distance", units.length)
"""
Distance between the centers of the particles.

Symbol:
    :code:`r`
"""

law = Eq(potential,
    4 * dispersion_energy * ((particle_size / distance)**12 - (particle_size / distance)**6))
r"""
:code:`U = 4 * e * ((sigma / r)^12 - (sigma / r)^6)`

Latex:
    .. math::
        U = 4 \varepsilon \left( \left( \frac{\sigma}{r} \right)^{12} - \left( \frac{\sigma}{r} \right)^6 \right)
"""

@validate_input(
    dispersion_energy_=dispersion_energy,
    particle_size_=particle_size,
    distance_=distance,
)
@validate_output(potential)
def calculate_potential(
    dispersion_energy_: Quantity,
    particle_size_: Quantity,
    distance_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        dispersion_energy: dispersion_energy_,
        particle_size: particle_size_,
        distance: distance_,
    })
    return Quantity(result)
