"""
Planet period squared is proportional to cube of semimajor axis
===============================================================

Also known as **Kepler's third law** of planetary motion, the law of periods relates
the period of rotation of any planet to the semi-major axis of its orbit.

**Links:**

#. `Physics LibreTexts. Kepler's Third Law, Derivation of Kepler's Third Law (5.6.23) <https://phys.libretexts.org/Bookshelves/University_Physics/Physics_(Boundless)/5%3A_Uniform_Circular_Motion_and_Gravitation/5.6%3A_Keplers_Laws>`__.
"""

from sympy import Eq, solve, pi, Symbol as SymSymbol
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    clone_as_symbol,
    symbols,
    quantities,
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

rotation_period = clone_as_symbol(symbols.period, positive=True)
"""
The planet's :symbols:`period` of rotation.
"""

attracting_mass = clone_as_symbol(symbols.mass, positive=True)
"""
:symbols:`mass` of the attracting body, e.g. the Sun in case of the solar system.
"""

semimajor_axis = clone_as_symbol(symbols.semimajor_axis, positive=True)
"""
:symbols:`semimajor_axis` of the planet's orbit. It is equal to the :symbols:`radius` in case of
a round orbit.
"""

law = Eq(
    rotation_period**2,
    4 * pi**2 / (quantities.gravitational_constant * attracting_mass) * semimajor_axis**3,
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive law from Newton's second law of motion

_angular_speed = SymSymbol("angular_speed", positive=True)
_radius = SymSymbol("radius", nonnegative=True)
_attracted_mass = newtons_second_law.mass

_acceleration_expr = centripetal_law.law.rhs.subs({
    centripetal_law.angular_speed: _angular_speed,
    centripetal_law.radius_of_curvature: _radius,
})

_force_expr = gravity_law.law.rhs.subs({
    gravity_law.distance_between_mass_centers: _radius,
    gravity_law.first_mass: attracting_mass,
    gravity_law.second_mass: _attracted_mass,
})

_newtons_eqn = newtons_second_law.law.subs({
    newtons_second_law.acceleration: _acceleration_expr,
    newtons_second_law.force: _force_expr,
    newtons_second_law.mass: _attracted_mass,
})

_angular_speed_expr = solve(_newtons_eqn, _angular_speed)[1]

# If the eccentricity of the planet's orbit is not too big, which is the case for most
# planets of the solar system, we can substitute the radius with the semi-major axis of
# the orbit.
_period_derived = period_law.law.rhs.subs(period_law.angular_frequency, _angular_speed_expr)
_period_derived = _period_derived.subs(_radius, semimajor_axis)

_period_from_law = solve(law, rotation_period)[0]

assert expr_equals(_period_derived, _period_from_law)


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
