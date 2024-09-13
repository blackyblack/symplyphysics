"""
Potential energy of radial planetary motion
===========================================

The total mechanical energy of the planet can be viewed as the sum of the kinetic and potential energy.
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

total_potential_energy = Symbol("total_potential_energy", units.energy)
r"""
The total potential energy of the planet.

Symbol:
    :code:`U_tot`

Latex:
    :math:`U_\text{tot}`
"""

gravitational_potential_energy = Symbol("gravitational_potential_energy", units.energy)
r"""
The potential energy of the planet due to the gravitational interaction of the planet and the star.

Symbol:
    :code:`U_gr`

Latex:
    :math:`U_\text{gr}`
"""

angular_momentum = Symbol("angular_momentum", units.length * units.momentum)
"""
The angular momentum of the planet.

Symbol:
    :code:`L`
"""

planetary_mass = symbols.mass
"""
The :attr:`~symplyphysics.symbols.mass` of the planet.
"""

distance = Symbol("distance", units.length)
"""
The distance between the star and the planet.

Symbol:
    :code:`r`
"""

law = Eq(
    total_potential_energy,
    gravitational_potential_energy + angular_momentum**2 / (2 * planetary_mass * distance**2),
)
r"""
:code:`U_tot = U_gr + L^2 / (2 * m * r^2)`

Latex:
    .. math::
        U_\text{tot} = U_\text{gr} + \frac{L^2}{2 m r^2}
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
