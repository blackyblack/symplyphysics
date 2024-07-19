"""
Potential energy of radial planetary motion
===========================================

The total energy of the planet can be viewed as the sum of the kinetic and potential energy.
The potential energy is in turn the sum of the potential energy due to the gravitational interaction
between the planet and the Sun, and the the energy of the tangential motion, which depends on the
planet's angular momentum.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
)

planetary_potential_energy = Symbol("planetary_potential_energy", units.energy)
"""
The total potential energy of the planet.

Symbol:
    V
"""

gravitational_potential_energy = Symbol("gravitational_potential_energy", units.energy)
"""
The potential energy of the planet due to the gravitational interaction of the planet and the star.

Symbol:
    U
"""

planetary_angular_momentum = Symbol("planetary_angular_momentum", units.length * units.momentum)
"""
The angular momentum of the planet.

Symbol:
    L
"""

planetary_mass = clone_symbol(symbols.basic.mass, "planetary_mass")
"""
The :attr:`symplyphysics.symbols.basic.mass` of the planet.

Symbol:
    m
"""

planetary_distance = Symbol("planetary_distance", units.length)
"""
The distance between the star and the planet.

Symbol:
    r
"""

law = Eq(
    planetary_potential_energy,
    gravitational_potential_energy + planetary_angular_momentum**2 / (2 * planetary_mass * planetary_distance**2),
)
r"""
V = U + L^2 / (2 * m * r^2)

Latex:
    :math:`V = U + \frac{L^2}{2 m r^2}`
"""


@validate_input(
    gravitational_potential_energy_=gravitational_potential_energy,
    planetary_angular_momentum_=planetary_angular_momentum,
    planetary_mass_=planetary_mass,
    planetary_distance_=planetary_distance,
)
@validate_output(planetary_potential_energy)
def calculate_planetary_potential_energy(
    gravitational_potential_energy_: Quantity,
    planetary_angular_momentum_: Quantity,
    planetary_mass_: Quantity,
    planetary_distance_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        gravitational_potential_energy: gravitational_potential_energy_,
        planetary_angular_momentum: planetary_angular_momentum_,
        planetary_mass: planetary_mass_,
        planetary_distance: planetary_distance_,
    })
    return Quantity(result)
