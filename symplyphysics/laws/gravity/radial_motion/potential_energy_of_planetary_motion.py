"""
Potential energy of radial planetary motion
===========================================

The total mechanical energy of the planet can be viewed as the sum of the kinetic and potential energy.
The potential energy is in turn the sum of the potential energy due to the gravitational interaction
between the planet and the Sun, and the the energy of the tangential motion, which depends on the
planet's angular momentum.

**Links:**

#. Sivukhin, D.V. (1979). __Obshchiy kurs fiziki__ [General course of Physics], vol. 1, p. 315.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

total_potential_energy = clone_as_symbol(
    symbols.potential_energy,
    display_symbol="U_tot",
    display_latex="U_\\text{tot}",
)
"""
The total :symbols:`potential_energy` of the planet.
"""

gravitational_potential_energy = clone_as_symbol(
    symbols.potential_energy,
    display_symbol="U_gr",
    display_latex="U_\\text{gr}",
)
"""
The :symbols:`potential_energy` of the planet due to the gravitational interaction of the planet and the star.
"""

angular_momentum = symbols.angular_momentum
"""
The :symbols:`angular_momentum` of the planet.
"""

planetary_mass = symbols.mass
"""
The :symbols:`mass` of the planet.
"""

distance = symbols.euclidean_distance
"""
The :symbols:`euclidean_distance` between the star and the planet.
"""

law = Eq(
    total_potential_energy,
    gravitational_potential_energy + angular_momentum**2 / (2 * planetary_mass * distance**2),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    gravitational_potential_energy_=gravitational_potential_energy,
    angular_momentum_=angular_momentum,
    planetary_mass_=planetary_mass,
    distance_=distance,
)
@validate_output(total_potential_energy)
def calculate_planetary_potential_energy(
    gravitational_potential_energy_: Quantity,
    angular_momentum_: Quantity,
    planetary_mass_: Quantity,
    distance_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        gravitational_potential_energy: gravitational_potential_energy_,
        angular_momentum: angular_momentum_,
        planetary_mass: planetary_mass_,
        distance: distance_,
    })
    return Quantity(result)
