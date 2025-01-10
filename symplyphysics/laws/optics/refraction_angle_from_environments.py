"""
Refraction angle from enviroments
=================================

If ray of light comes from one medium to another, it refracts.
Incoming ray, refracted ray and perpendicular to medium boundary are within the same plane.

**Conditions:**

#. Light is monochromatic, as the refractive index depends on the light frequency.
#. The refracting medium is uniform, so that the refracting index does not change over the ray path.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Snell%27s_law>`__.

..
    TODO rename files
"""

from sympy import (Eq, solve, diff, sin, pi, sqrt)
from symplyphysics import (units, Quantity, angle_type,
    validate_input, validate_output, symbols, clone_as_symbol)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import refractive_index_is_wave_speeds_ratio as refractive_index_definition
from symplyphysics.laws.kinematics import position_via_constant_speed_and_time as distance_law
from symplyphysics.laws.geometry import planar_projection_is_cosine as projection_law
from symplyphysics.core.symbols.quantities import scale_factor

incidence_refractive_index = clone_as_symbol(symbols.relative_refractive_index, subscript="1")
"""
:symbols:`relative_refractive_index` of the medium in which the indicent ray travels.
"""

resulting_refractive_index = clone_as_symbol(symbols.relative_refractive_index, subscript="2")
"""
:symbols:`relative_refractive_index` of the medium in which the refracted ray travels.
"""

incidence_angle = clone_as_symbol(symbols.angle, subscript="1")
"""
:symbols:`angle` of incidence.
"""

refraction_angle = clone_as_symbol(symbols.angle, subscript="2")
"""
:symbols:`angle` of refraction.
"""

law = Eq(incidence_refractive_index * sin(incidence_angle),
    resulting_refractive_index * sin(refraction_angle))
"""
:laws:symbol::

:laws:latex::
"""

# Derive the same Snell's law from Fermat's principle
# TODO: prefix variables used in proof with underscore

# Fermat's principle states that the path taken by a ray between two given points is the path that can be traveled in the least time.
# For derivation, we will use the notation in figure https://en.wikipedia.org/wiki/Fermat%27s_principle#/media/File:Fermat_Snellius.svg
# First, let's find an expression for the time it takes light to travel from A to B in terms of _x - point of refraction at the boundary of the medium.
_travel_time = symbols.time
_outer_travel_distance = clone_as_symbol(symbols.euclidean_distance, subscript="1")
_medium_travel_distance = clone_as_symbol(symbols.euclidean_distance, subscript="2")
_outer_velocity = clone_as_symbol(symbols.speed, subscript="1")
_medium_velocity = clone_as_symbol(symbols.speed, subscript="2")

_distance_law_outer_eq = distance_law.law.subs({
    distance_law.final_position: _outer_travel_distance,
    distance_law.speed: _outer_velocity,
    distance_law.initial_position: 0
})
_distance_law_medium_eq = distance_law.law.subs({
    distance_law.final_position: _medium_travel_distance,
    distance_law.speed: _medium_velocity,
    distance_law.initial_position: 0
})
_outer_travel_time = solve(_distance_law_outer_eq, distance_law.time)[0]
_medium_travel_time = solve(_distance_law_medium_eq, distance_law.time)[0]

_travel_time = _outer_travel_time + _medium_travel_time

# Lengths _outer_travel_distance (l1) and _medium_travel_distance (l2) can be expressed in terms of _x.
_x = symbols.position
_a = clone_as_symbol(symbols.length, display_symbol="a")
_b = clone_as_symbol(symbols.length, display_symbol="b")
_d = clone_as_symbol(symbols.length, display_symbol="d")
_outer_travel_distance_cartesian_eq = Eq(_outer_travel_distance, sqrt(_x**2 + _a**2))
_medium_travel_distance_cartesian_eq = Eq(_medium_travel_distance, sqrt((_d - _x)**2 + _b**2))

# It is also possible to express _outer_travel_distance and _medium_travel_distance in terms of the angles of incidence (alpha) and refraction (beta).

# Use (pi / 2 - angle) to obtain vertical projection instead of horizontal
_projection_incidence_eq = projection_law.law.subs({
    projection_law.vector_angle: pi / 2 - incidence_angle,
    projection_law.vector_length: _outer_travel_distance,
    projection_law.projection: _x
})
_projection_incidence_distance = solve(_projection_incidence_eq, _outer_travel_distance)[0]
_outer_travel_distance_polar_eq = Eq(_outer_travel_distance, _projection_incidence_distance)
_projection_refraction_eq = projection_law.law.subs({
    projection_law.vector_angle: pi / 2 - refraction_angle,
    projection_law.vector_length: _medium_travel_distance,
    projection_law.projection: _d - _x
})
_projection_refraction_distance = solve(_projection_refraction_eq, _medium_travel_distance)[0]
_medium_travel_distance_polar_eq = Eq(_medium_travel_distance, _projection_refraction_distance)

# NOTE: derivative of the path equals to zero does not always correspond to minumum travel time. It
#       can also be a maximum time. Here we prove that Snell's law is applicable to cases
#       with extreme travel time.
# NOTE: light travel path is known to have minimum time, which is an extreme time case.

# Substitute l1, l2 in _travel_time, differentiate by _x and equate the derivative to zero - extreme time case.
_travel_time_on_cartesian = _travel_time.subs({
    _outer_travel_distance_cartesian_eq.lhs: _outer_travel_distance_cartesian_eq.rhs,
    _medium_travel_distance_cartesian_eq.lhs: _medium_travel_distance_cartesian_eq.rhs
})
_extreme_time_case = Eq(diff(_travel_time_on_cartesian, _x), 0)

# Let's get the same expression using the sines of the angles.
_extreme_time_case = _extreme_time_case.subs({
    _outer_travel_distance_cartesian_eq.rhs: _outer_travel_distance_cartesian_eq.lhs,
    _medium_travel_distance_cartesian_eq.rhs: _medium_travel_distance_cartesian_eq.lhs
})
_extreme_time_case = _extreme_time_case.subs({
    _outer_travel_distance_polar_eq.lhs: _outer_travel_distance_polar_eq.rhs,
    _medium_travel_distance_polar_eq.lhs: _medium_travel_distance_polar_eq.rhs
}).simplify()

# Finally, let's use the definition of the refractive index as the ratio of the speed of light in the medium to that in _a reference medium (vacuum).
_outer_refraction_definition = refractive_index_definition.definition.subs(
    {refractive_index_definition.relative_refractive_index: incidence_refractive_index})
_medium_refreaction_definition = refractive_index_definition.definition.subs(
    {refractive_index_definition.relative_refractive_index: resulting_refractive_index})
_outer_refraction_velocity = solve(_outer_refraction_definition,
    refractive_index_definition.refracted_wave_speed,
    dict=True)[0][refractive_index_definition.refracted_wave_speed]
_medium_refraction_velocity = solve(_medium_refreaction_definition,
    refractive_index_definition.refracted_wave_speed,
    dict=True)[0][refractive_index_definition.refracted_wave_speed]

_extreme_time_case = _extreme_time_case.subs({
    _outer_velocity: _outer_refraction_velocity,
    _medium_velocity: _medium_refraction_velocity
})
_resulting_expression = Eq(
    incidence_refractive_index * sin(incidence_angle),
    solve(_extreme_time_case, incidence_refractive_index * sin(incidence_angle),
    dict=True)[0][incidence_refractive_index * sin(incidence_angle)])

# Verify that the resulting_expression corresponds to the law.
assert expr_equals(law.lhs, _resulting_expression.lhs)
assert expr_equals(law.rhs, _resulting_expression.rhs)


@validate_input(incidence_angle_=incidence_angle,
    incidence_refractive_index_=incidence_refractive_index,
    resulting_refractive_index_=resulting_refractive_index)
@validate_output(refraction_angle)
def calculate_refraction_angle(incidence_angle_: Quantity | float,
    incidence_refractive_index_: float, resulting_refractive_index_: float) -> Quantity:
    #HACK: sympy angles are always in radians
    incidence_angle_radians = scale_factor(incidence_angle_)
    # Check for boundary conditions
    assert incidence_angle_radians <= pi / 2
    assert incidence_angle_radians >= -pi / 2
    solutions = solve(law, refraction_angle, dict=True)
    result_expr = solutions[0][refraction_angle]
    angle_applied = result_expr.subs({
        incidence_angle: incidence_angle_radians,
        incidence_refractive_index: incidence_refractive_index_,
        resulting_refractive_index: resulting_refractive_index_
    })
    #HACK: there are 2 solutions for refraction_angle: pi - asin() and asin(). We choose former and switch to latter
    #      if resulting angle does not pass boundary conditions.
    if (angle_applied > pi / 2 or angle_applied < -pi / 2):
        result_expr = solutions[1][refraction_angle]
        angle_applied = result_expr.subs({
            incidence_angle: incidence_angle_radians,
            incidence_refractive_index: incidence_refractive_index_,
            resulting_refractive_index: resulting_refractive_index_
        })

    # Check for boundary conditions
    assert angle_applied <= pi / 2
    assert angle_applied >= -pi / 2
    #HACK: angle type is automatically detected as dimensionless. Force it to angle.
    return Quantity(angle_applied * units.radian, dimension=angle_type)
