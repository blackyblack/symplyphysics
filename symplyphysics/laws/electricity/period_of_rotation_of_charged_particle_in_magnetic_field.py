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
"""

from sympy import (Eq, solve, pi, symbols as sympy_symbols)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, symbols, Vector, vector_magnitude)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity.vector import lorentz_force_via_electromagnetic_field as _lorentz_law
from symplyphysics.laws.dynamics import acceleration_is_force_over_mass as _newtons_second_law
from symplyphysics.laws.kinematics import (
    centripetal_acceleration_via_angular_speed_and_radius as _centripetal_acceleration_law,
    speed_via_angular_speed_and_radius as _speed_law,
)
from symplyphysics.definitions import period_from_angular_frequency as _period_law

period = Symbol("period", units.time)
"""
Period of the particle's rotation.

Symbol:
    :code:`T`
"""

mass = symbols.basic.mass
"""
:attr:`~symplyphysics.symbols.basic.mass` of the particle.
"""

charge = Symbol("charge", units.charge)
"""
Charge of the particle.

Symbol:
    :code:`q`
"""

magnetic_flux_density = Symbol("magnetic_flux_density", units.magnetic_density)
"""
Magnitude of magnetic flux density.

Symbol:
    :code:`B`
"""

law = Eq(period, 2 * pi * mass / (charge * magnetic_flux_density))
r"""
:code:`T = 2 * pi * m / (q * B)`

Latex:
    .. math::
        T = 2 \pi \frac{m}{q B}
"""

# Derive law from Lorentz force

_charge = sympy_symbols("charge", positive=True)
_angular_frequency = sympy_symbols("angular_frequency", positive=True)
_radius = sympy_symbols("radius", positive=True)

_speed_expr = _speed_law.law.rhs.subs({
    _speed_law.angular_speed: _angular_frequency,
    _speed_law.radius_of_curvature: _radius,
})

_velocity_vec = Vector([_speed_expr, 0, 0])

_magnetic_field = sympy_symbols("magnetic_field", positive=True)
_magnetic_field_vec = Vector([0, _magnetic_field, 0])

_force_vec = _lorentz_law.lorentz_force_law(
    electric_field_=Vector([0, 0, 0]),
    magnetic_flux_density_=_magnetic_field_vec,
    velocity_=_velocity_vec,
).subs(
    _lorentz_law.charge, _charge
)

_force_expr = vector_magnitude(_force_vec)

_acceleration_via_force = _newtons_second_law.law.rhs.subs({
    _newtons_second_law.force: _force_expr,
    _newtons_second_law.mass: mass,
})

_acceleration_via_frequency = _centripetal_acceleration_law.law.rhs.subs({
    _centripetal_acceleration_law.angular_speed: _angular_frequency,
    _centripetal_acceleration_law.radius_of_curvature: _radius,
})

_angular_frequency_expr = solve(
    eq := Eq(_acceleration_via_frequency, _acceleration_via_force),
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
