"""
Potential energy of radial planetary motion
===========================================

TODO
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
TODO

Symbol:
    V
"""

gravitational_potential_energy = Symbol("gravitational_potential_energy", units.energy)
"""
TODO

:doc:`../gravitational_potential_energy.py`

Symbol:
    U
"""

planetary_angular_momentum = Symbol("planetary_angular_momentum", units.length * units.momentum)
"""
TODO

Symbol:
    L
"""

planetary_mass = clone_symbol(symbols.basic.mass, "planetary_mass")
"""
TODO

Symbol:
    m
"""

planetary_distance = Symbol("planetary_distance", units.length)
"""
TODO

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
