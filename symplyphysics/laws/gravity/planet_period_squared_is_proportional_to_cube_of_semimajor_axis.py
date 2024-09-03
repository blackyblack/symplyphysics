from sympy import Eq, solve, pi, Symbol as SymSymbol
from sympy.physics.units import gravitational_constant
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    clone_symbol,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import (
    period_from_angular_frequency as period_law,)
from symplyphysics.laws.dynamics import (
    acceleration_is_force_over_mass as newtons_second_law,)
from symplyphysics.laws.gravity import (
    gravity_force_from_mass_and_distance as gravity_law,)
from symplyphysics.laws.kinematics import (
    centripetal_acceleration_via_angular_speed_and_radius as centripetal_law,)

# Description
## Also known as Kepler's third law of planetary motion, the law of periods relates
## the period of rotation of any planet to the semi-major axis of its orbit.

# Law: T**2 = 4*pi**2 / (G*M) * a**3
## T - planet's period of rotation
## G - gravitational constant
## M - mass of the attracting body, the Sun in case of the solar system
## a - semi-major axis of the planet's orbit, its radius in case of a round orbit

rotation_period = Symbol("rotation_period", units.time, positive=True)
attracting_mass = clone_symbol(symbols.basic.mass, positive=True)
semimajor_axis = Symbol("semimajor_axis", units.length, positive=True)

law = Eq(
    rotation_period**2,
    4 * pi**2 / (gravitational_constant * attracting_mass) * semimajor_axis**3,
)

# Derive law from Newton's second law of motion

_angular_speed = SymSymbol("angular_speed", positive=True)
_radius = SymSymbol("radius", nonnegative=True)
_attracted_mass = newtons_second_law.mass

_acceleration_expr = centripetal_law.law.rhs.subs({
    centripetal_law.angular_speed: _angular_speed,
    centripetal_law.radius_of_curvature: _radius,
})

_force_expr = gravity_law.law.rhs.subs({
    gravity_law.first_mass: attracting_mass,
    gravity_law.distance_between_mass_centers: _radius,
    gravity_law.second_mass: _attracted_mass,
})

_newtons_eqn = newtons_second_law.law.subs({
    newtons_second_law.acceleration: _acceleration_expr,
    newtons_second_law.force: _force_expr,
})

_angular_speed_expr = solve(_newtons_eqn, _angular_speed)[1]

# If the eccentricity of the planet's orbit is not too big, which is the case for most
# planets of the solar system, we can substitute the radius with the semi-major axis of
# the orbit.
_period_derived = period_law.law.rhs.subs({
    period_law.angular_frequency: _angular_speed_expr,
    _radius: semimajor_axis,
})

_period_from_law = solve(law, rotation_period)[0]

assert expr_equals(_period_derived, _period_from_law)


def print_law() -> str:
    return print_expression(law)


@validate_input(attracting_mass_=attracting_mass, semimajor_axis_=semimajor_axis)
@validate_output(rotation_period)
def calculate_rotation_period(
    attracting_mass_: Quantity,
    semimajor_axis_: Quantity,
) -> Quantity:
    result = solve(law, rotation_period)[0].subs({
        attracting_mass: attracting_mass_,
        semimajor_axis: semimajor_axis_,
    })
    return Quantity(result)
