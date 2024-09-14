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
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
)

planetary_mass = symbols.mass
"""
The :attr:`~symplyphysics.symbols.mass` of the planet.
"""

radial_speed = Symbol("radial_speed", units.velocity)
"""
The projection of the velocity vector in the radial direction.

Symbol:
    :code:`v_r`

Latex:
    :math:`v_r`
"""

potential_energy = Symbol("potential_energy", units.energy)
"""
The :doc:`potential energy <laws.gravity.radial_motion.potential_energy_of_planetary_motion>` of the planet.

Symbol:
    :code:`U`
"""

total_mechanical_energy = Symbol("total_mechanical_energy", units.energy)
"""Total mechanical energy of the planet, assumed to be constant.

Symbol:
    :code:`E`
"""

law = Eq(planetary_mass * radial_speed**2 / 2 + potential_energy, total_mechanical_energy)
r"""
:code:`m * v_r**2 / 2 + U = E`

Latex:
    .. math::
        \frac{1}{2} m v_r^2 + U = E
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
