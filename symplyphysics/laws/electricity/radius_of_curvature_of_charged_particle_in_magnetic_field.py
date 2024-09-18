"""
Radius of curvature of charged particle in magnetic field
=========================================================

When a charged particle enters a magnetic field, it experiences an :doc:`electromagnetic
force <laws.electricity.vector.lorentz_force_via_electromagnetic_field>` upon itself.
In the absence of the electric field, the particle starts moving in a circular orbit. The
radius of curvature of the particle's orbit is determined by the mass, speed, and charge
of the particle as well as by the magnetic flux density.

**Conditions:**

#. The particle's speed and the magnetic field are perpendicular to each other.
#. The magnetic field is uniform.
#. The electric field is zero.
"""

from sympy import (Eq, solve, pi)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity import period_of_rotation_of_charged_particle_in_magnetic_field as period_law
from symplyphysics.laws.kinematics import position_via_constant_speed_and_time as distance_law

radius_of_curvature = symbols.radius_of_curvature
"""
Radius of curvature of the particle's orbit.
"""

mass = symbols.mass
"""
:symbols:`mass` of the particle.
"""

speed = symbols.speed
"""
Speed of the particle.
"""

charge = symbols.charge
"""
Charge of the particle.
"""

magnetic_flux_density = symbols.magnetic_flux_density
"""
Magnitude of magnetic flux density.
"""

law = Eq(radius_of_curvature, mass * speed / (charge * magnetic_flux_density))
"""
:laws:symbol::

:laws:latex::
"""

# This law might be derived via period of a charged particle in a magnetic field and distance from constant speed.

_distance_law_applied = distance_law.law.subs({
    distance_law.initial_position: 0,
    distance_law.final_position: 2 * pi * radius_of_curvature,
    distance_law.speed: speed,
})
# Period is a time taken to cover a full circle.
_period_derived = solve(_distance_law_applied, distance_law.time,
    dict=True)[0][distance_law.time]

_law_applied = period_law.law.subs({
    period_law.mass: mass,
    period_law.charge: charge,
    period_law.magnetic_flux_density: magnetic_flux_density,
    period_law.period: _period_derived,
})
_radius_derived = solve(_law_applied, radius_of_curvature, dict=True)[0][radius_of_curvature]

# Check if derived radius_of_curvature is same as declared.
assert expr_equals(_radius_derived, law.rhs)


@validate_input(mass_=mass, velocity_=speed, induction_=magnetic_flux_density, charge_=charge)
@validate_output(radius_of_curvature)
def calculate_radius(mass_: Quantity, velocity_: Quantity, induction_: Quantity,
    charge_: Quantity) -> Quantity:
    result_expr = solve(law, radius_of_curvature, dict=True)[0][radius_of_curvature]
    result_expr = result_expr.subs({
        mass: mass_,
        speed: velocity_,
        magnetic_flux_density: induction_,
        charge: charge_,
    })
    return Quantity(result_expr)
