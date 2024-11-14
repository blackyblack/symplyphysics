"""
Gravitational potential energy
==============================

The **gravitational potential energy** of a system of two particles is the negative of the work
that would be done by the gravitational force of either particle acting on the other if the
particles were brought together from infinity to the given distance.

**Notation:**

#. :quantity_notation:`gravitational_constant`.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    clone_as_symbol,
    symbols,
    quantities,
)

gravitational_potential_energy = symbols.potential_energy
"""
Gravitational :symbols:`potential_energy`.
"""

first_mass = clone_as_symbol(symbols.mass, subscript="1")
"""
:symbols:`mass` of the first particle.
"""

second_mass = clone_as_symbol(symbols.mass, subscript="2")
"""
:symbols:`mass` of the second particle.
"""

distance_between_mass_centers = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` between mass centers of the particles.
"""


law = Eq(gravitational_potential_energy,
    -1 * quantities.gravitational_constant * first_mass * second_mass / distance_between_mass_centers)


@validate_input(
    first_mass_=first_mass,
    second_mass_=second_mass,
    distance_between_mass_centers_=distance_between_mass_centers,
)
@validate_output(gravitational_potential_energy)
def calculate_gravitational_potential_energy(
    first_mass_: Quantity,
    second_mass_: Quantity,
    distance_between_mass_centers_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        first_mass: first_mass_,
        second_mass: second_mass_,
        distance_between_mass_centers: distance_between_mass_centers_,
    })
    return Quantity(result)
