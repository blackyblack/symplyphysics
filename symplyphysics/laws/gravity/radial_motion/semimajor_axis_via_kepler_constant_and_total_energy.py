r"""
Semimajor axis via Kepler's constant and total energy
=====================================================

The semi-major axis of an orbiting planet depends on the Kepler's constant of
the star---planet system and the total energy of the planet per unit of its mass.

**Notes:**

#. Works for both elliptical (:math:`\varepsilon < 0`) and hyperbolical
   (:math:`\varepsilon > 0`) orbits.
"""

from sympy import Eq, pi
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

semimajor_axis = Symbol("semimajor_axis", units.length)
"""
The semi-major axis of the planet's orbit.

Symbol:
    a
"""

kepler_constant = Symbol("kepler_constant", units.length**3 / units.time**2)
r"""
The Kepler's constant, whose value is determined by the mass of the orbited star.

Symbol:
    K

Latex:
    :math:`\mathfrak{K}`
"""

total_energy_per_unit_mass = Symbol("total_energy_per_unit_mass", units.energy / units.mass)
r"""
The total energy of the planet per unit of its :attr:`~symplyphysics.symbols.basic.mass`.
Can be negative or positive.

Symbol:
    epsilon

Latex:
    :math:`\varepsilon`
"""

law = Eq(semimajor_axis, 2 * pi**2 * kepler_constant / abs(total_energy_per_unit_mass))
r"""
a = 2 * pi^2 * K / abs(epsilon)

Latex:
    .. math::
        a = \frac{2 \pi^2 \mathfrak{K}}{|\varepsilon|}
"""


@validate_input(
    kepler_constant_=kepler_constant,
    total_energy_per_unit_mass_=total_energy_per_unit_mass,
)
@validate_output(semimajor_axis)
def calculate_semimajor_axis(
    kepler_constant_: Quantity,
    total_energy_per_unit_mass_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        kepler_constant: kepler_constant_,
        total_energy_per_unit_mass: total_energy_per_unit_mass_,
    })
    return Quantity(result)
