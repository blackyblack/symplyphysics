"""
Second escape velocity
======================

The second cosmic velocity is the lowest velocity that must be given to an object starting from the surface of
a celestial body, the mass of which is negligible compared to the mass of a celestial body, in order to overcome
the gravitational attraction of this celestial body and leave a closed orbit around it.

**Notation:**

#. :quantity_notation:`gravitational_constant`.
"""

from sympy import Eq, solve, sqrt
from symplyphysics import symbols, Quantity, validate_input, validate_output, quantities
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.gravity import gravitational_potential_energy as potential_energy_law
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_speed as kinetic_energy_law
from symplyphysics.laws.conservation import mechanical_energy_after_equals_to_mechanical_energy_before as conservation_law

speed = symbols.speed
"""
Second escape :symbols:`speed`.
"""

planet_mass = symbols.mass
"""
Planetary :symbols:`mass`.
"""

planet_radius = symbols.radius
"""
Planetary :symbols:`radius`.
"""

height = symbols.height
"""
Elevation (see :symbols:`height`) of the body above the planet's surface.
"""

law = Eq(speed, sqrt(2 * quantities.gravitational_constant * planet_mass / (planet_radius + height)))
"""
:laws:symbol::

:laws:latex::
"""

# This law might be derived via "gravitational_potential_energy" law, "kinetic_energy_from_mass_and_velocity" law,
# "mechanical_energy_after_equals_to_mechanical_energy_before" law.

# The radius of the orbit consists of the radius of the planet and the height above its surface.
_potential_energy_law_applied = potential_energy_law.law.subs({
    potential_energy_law.first_mass: planet_mass,
    potential_energy_law.distance_between_mass_centers: planet_radius + height,
})
_potential_energy_derived = solve(_potential_energy_law_applied,
    potential_energy_law.gravitational_potential_energy,
    dict=True)[0][potential_energy_law.gravitational_potential_energy]

_kinetic_energy_law_applied = kinetic_energy_law.law.subs({
    kinetic_energy_law.mass: potential_energy_law.second_mass,
    kinetic_energy_law.speed: speed,
})
_kinetic_energy_derived = solve(_kinetic_energy_law_applied,
    kinetic_energy_law.kinetic_energy,
    dict=True)[0][kinetic_energy_law.kinetic_energy]

# Here we consider the fall of a body onto a planet from infinity. In this case, the law of conservation of energy will
# look like: kinetic energy minus potential. But in the law "gravitational_potential_energy", the potential energy has
# already been introduced with a minus. And in the law "mechanical_energy_after_equals_to_mechanical_energy_before",
# the conservation law looks like: E1=E2. Therefore, a minus sign is placed before the potential energy.
_conservation_law_applied = conservation_law.law.subs({
    conservation_law.mechanical_energy(conservation_law.time_after): _kinetic_energy_derived,
    conservation_law.mechanical_energy(conservation_law.time_before): -1 * _potential_energy_derived,
})
# The proof considers the fall of a body from infinity to the earth. But at the same time, in most cases, I consider
# this speed as directed from the planet from which the body starts, which corresponds to the plus sign. Therefore,
# the first solution is ignored.
_velocity_derived = solve(_conservation_law_applied, speed, dict=True)[1][speed]

# Check if derived velocity is same as declared.
assert expr_equals(_velocity_derived, law.rhs)


@validate_input(planet_mass_=planet_mass, radius_=planet_radius, height_=height)
@validate_output(speed)
def calculate_velocity(planet_mass_: Quantity, radius_: Quantity, height_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, speed, dict=True)[0][speed]
    result_expr = result_velocity_expr.subs({
        planet_mass: planet_mass_,
        planet_radius: radius_,
        height: height_
    })
    return Quantity(result_expr)
