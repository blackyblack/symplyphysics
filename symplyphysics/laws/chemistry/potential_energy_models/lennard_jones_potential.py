"""
Lennard-Jones potential
=======================

The Lennard-Jones (LJ) potential is an intermolecular pair potential. It models soft repulsive and
attractive (van der Waals) interactions, and it describes electronically neutral atoms or molecules.
It is a simple yet realistic model to describe intermolecular interactions: two particles repel each
other at a very close distance, attract each other at moderate distance, and do not interact at infinite
distance.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Lennard-Jones_potential#Overview>`__.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

potential = symbols.potential_energy
"""
:symbols:`potential_energy` of the configuration.
"""

dispersion_energy = clone_as_symbol(symbols.energy, display_symbol="e", display_latex="\\varepsilon")
"""
Depth of the potential well, also referred to as *dispersion energy*. See :symbols:`energy`.
"""

particle_size = clone_as_symbol(symbols.radius, display_symbol="sigma", display_latex="\\sigma")
"""
Distance at which potential is zero, also referred to as *particle size*.
"""

distance = clone_as_symbol(symbols.euclidean_distance, display_symbol="r", display_latex="r")
"""
:symbols:`euclidean_distance` between the centers of the particles.
"""

law = Eq(potential,
    4 * dispersion_energy * ((particle_size / distance)**12 - (particle_size / distance)**6))
"""
:laws:symbol::

:laws:latex::
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
