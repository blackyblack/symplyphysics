"""
Period of rotation of charged particle in magnetic field
========================================================

When a charged particle enters a magnetic field, it experiences an :doc:`electromagnetic
force <laws.electricity.vector.lorentz_force_via_electromagnetic_field>` upon itself.
In the absence of the electric field, the particle starts moving in a circular orbit. The
period of the particle's rotation is determined by the mass and charge of the particle
as well as by the magnetic field and it does not depend on the particle speed.

**Conditions:**

#. The particle's speed and the magnetic field are perpendicular to each other.
#. The magnetic field is uniform.
#. The electric field is zero.

**Links:**

#. `Physics LibreTexts, formula 11.4.2 <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_II_-_Thermodynamics_Electricity_and_Magnetism_(OpenStax)/11%3A_Magnetic_Forces_and_Fields/11.04%3A_Motion_of_a_Charged_Particle_in_a_Magnetic_Field>`__.
"""

from sympy import Eq, solve, pi, symbols as sympy_symbols
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity.vector import lorentz_force_via_electromagnetic_field as _lorentz_law
from symplyphysics.laws.dynamics import acceleration_is_force_over_mass as _newtons_second_law
from symplyphysics.laws.kinematics import (
    centripetal_acceleration_via_angular_speed_and_radius as _centripetal_acceleration_law,
    speed_via_angular_speed_and_radius as _speed_law,
)
from symplyphysics.definitions import period_from_angular_frequency as _period_law

from symplyphysics.core.experimental.vectors import VectorNorm
from symplyphysics.core.experimental.coordinate_systems import CoordinateVector, CARTESIAN

period = symbols.period
"""
:symbols:`period` of the particle's rotation.
"""

mass = symbols.mass
"""
:symbols:`mass` of the particle.
"""

charge = symbols.charge
"""
:symbols:`charge` of the particle.
"""

magnetic_flux_density = symbols.magnetic_flux_density
"""
Magnitude of :symbols:`magnetic_flux_density`.
"""

law = Eq(period, 2 * pi * mass / (charge * magnetic_flux_density))
"""
:laws:symbol::

:laws:latex::
"""

# Derive law from Lorentz force

_charge = sympy_symbols("charge", positive=True)
_angular_frequency = sympy_symbols("angular_frequency", positive=True)
_radius = sympy_symbols("radius", positive=True)

_speed_expr = _speed_law.law.rhs.subs({
    _speed_law.angular_speed: _angular_frequency,
    _speed_law.radius_of_curvature: _radius,
})

_velocity_vec_ = CoordinateVector([_speed_expr, 0, 0], CARTESIAN)

_magnetic_field = sympy_symbols("magnetic_field", positive=True)
_magnetic_field_vec_ = CoordinateVector([0, _magnetic_field, 0], CARTESIAN)

_force_vec_ = _lorentz_law.law.rhs.subs({
    _lorentz_law.charge: _charge,
    _lorentz_law.electric_field: 0,  # NOTE: see condition 3 of this law
    _lorentz_law.velocity: _velocity_vec_,
    _lorentz_law.magnetic_flux_density: _magnetic_field_vec_,
})

_force_expr = VectorNorm(_force_vec_).doit()

_acceleration_via_force = _newtons_second_law.law.rhs.subs({
    _newtons_second_law.force: _force_expr,
    _newtons_second_law.mass: mass,
})

_acceleration_via_frequency = _centripetal_acceleration_law.law.rhs.subs({
    _centripetal_acceleration_law.angular_speed: _angular_frequency,
    _centripetal_acceleration_law.radius_of_curvature: _radius,
})

_angular_frequency_expr = solve(
    Eq(_acceleration_via_frequency, _acceleration_via_force),
    _angular_frequency,
)[0]

_period_expr = _period_law.law.rhs.subs({
    _period_law.angular_frequency: _angular_frequency_expr
}).subs({
    _charge: charge,
    _magnetic_field: magnetic_flux_density,
})

assert expr_equals(_period_expr, law.rhs)


@validate_input(mass_=mass, charge_=charge, induction_=magnetic_flux_density)
@validate_output(period)
def calculate_period(mass_: Quantity, charge_: Quantity, induction_: Quantity) -> Quantity:
    result_period_expr = solve(law, period, dict=True)[0][period]
    result_expr = result_period_expr.subs({
        mass: mass_,
        charge: charge_,
        magnetic_flux_density: induction_
    })
    return Quantity(result_expr)
