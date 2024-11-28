"""
Radial kinetic energy plus potential energy is constant
=======================================================

The total energy of the planet is composed of its radial kinetic energy and potential
energy. Note that the sign of the total energy determines the type of the planet's orbit:

#. If :math:`E < 0`, the orbit is elliptical.
#. If :math:`E = 0`, the orbit is parabolical.
#. If :math:`E > 0`, the orbit is hyperbolical.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

planetary_mass = symbols.mass
"""
The :symbols:`mass` of the planet.
"""

radial_speed = clone_as_symbol(symbols.speed, subscript="r")
"""
The projection of the velocity vector in the radial direction. See :symbols:`speed`.
"""

potential_energy = symbols.potential_energy
"""
The :doc:`potential energy <laws.gravity.radial_motion.potential_energy_of_planetary_motion>` of the planet.
See :symbols:`potential_energy`.
"""

total_mechanical_energy = symbols.mechanical_energy
"""Total :symbols:`mechanical_energy` of the planet, assumed to be constant.
"""

law = Eq(planetary_mass * radial_speed**2 / 2 + potential_energy, total_mechanical_energy)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    planetary_mass_=planetary_mass,
    radial_speed_=radial_speed,
    potential_energy_=potential_energy,
)
@validate_output(total_mechanical_energy)
def calculate_total_energy(
    planetary_mass_: Quantity,
    radial_speed_: Quantity,
    potential_energy_: Quantity,
) -> Quantity:
    result = law.lhs.subs({
        planetary_mass: planetary_mass_,
        radial_speed: radial_speed_,
        potential_energy: potential_energy_,
    })
    return Quantity(result)
