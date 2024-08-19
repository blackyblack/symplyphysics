from sympy import (Eq, solve, diff, sin, pi, sqrt)
from symplyphysics import (units, Quantity, Symbol, print_expression, dimensionless, angle_type,
    validate_input, validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import refractive_index_is_wave_speeds_ratio as refractive_index_definition
from symplyphysics.laws.kinematics import distance_from_constant_velocity as distance_law
from symplyphysics.laws.geometry import planar_projection_is_cosine as projection_law
from symplyphysics.core.symbols.quantities import scale_factor

# Description
## If ray of light comes from one medium to another, it refracts.
## Incoming ray, refracted ray and perpendicular to medium boundary are within the same plane.

## Law: n1 * sin(alpha) = n2 * sin(beta)
## Where:
## n1 is refractive index of first medium.
## n2 is refractive index of second medium.
## alfa is incoming angle. Angle is measured between normal vector to a medium boundary plane directed outside of refracting medium and incoming ray of light.
## beta is refraction angle. Angle is measured between normal vector to a medium boundary plane directed inside refracting medium and refracted ray of light.

# Conditions
## - alfa is in [-pi/2, pi/2] range, meaning that incoming ray comes from the outer medium.
## - beta is in [-pi/2, pi/2] range, meaning that refracted ray stays in the refracting medium.
## - light is monochromic, as refactive index depends on the light frequency.
## - refracting medium is uniform, so refracting index does not change over ray path.

incidence_refractive_index = Symbol("incidence_refractive_index", dimensionless)
resulting_refractive_index = Symbol("resulting_refractive_index", dimensionless)
incidence_angle = Symbol("incidence_angle", angle_type)
refraction_angle = Symbol("refraction_angle", angle_type)

law = Eq(incidence_refractive_index * sin(incidence_angle),
    resulting_refractive_index * sin(refraction_angle))

# Derive the same Snell's law from Fermat's principle

# Fermat's principle states that the path taken by a ray between two given points is the path that can be traveled in the least time.
# For derivation, we will use the notation in figure https://en.wikipedia.org/wiki/Fermat%27s_principle#/media/File:Fermat_Snellius.svg
# First, let's find an expression for the time it takes light to travel from A to B in terms of x - point of refraction at the boundary of the medium.
travel_time = Symbol("travel_time", units.time)
outer_travel_distance = Symbol("outer_travel_distance", units.length)
medium_travel_distance = Symbol("medium_travel_distance", units.length)
outer_velocity = Symbol("outer_velocity", units.velocity)
medium_velocity = Symbol("medium_velocity", units.velocity)

distance_law_outer_eq = distance_law.law.subs({
    distance_law.distance(distance_law.movement_time): outer_travel_distance,
    distance_law.constant_velocity: outer_velocity,
    distance_law.initial_position: 0
})
distance_law_medium_eq = distance_law.law.subs({
    distance_law.distance(distance_law.movement_time): medium_travel_distance,
    distance_law.constant_velocity: medium_velocity,
    distance_law.initial_position: 0
})
outer_travel_time = solve(distance_law_outer_eq, distance_law.movement_time)[0]
medium_travel_time = solve(distance_law_medium_eq, distance_law.movement_time)[0]

travel_time = outer_travel_time + medium_travel_time

# Lengths outer_travel_distance (l1) and medium_travel_distance (l2) can be expressed in terms of x.
x = Symbol("x", units.length)
a = Symbol("a", units.length)
b = Symbol("b", units.length)
d = Symbol("d", units.length)
outer_travel_distance_cartesian_eq = Eq(outer_travel_distance, sqrt(x**2 + a**2))
medium_travel_distance_cartesian_eq = Eq(medium_travel_distance, sqrt((d - x)**2 + b**2))

# It is also possible to express outer_travel_distance and medium_travel_distance in terms of the angles of incidence (alpha) and refraction (beta).

# Use (pi / 2 - angle) to obtain vertical projection instead of horizontal
projection_incidence_eq = projection_law.law.subs({
    projection_law.vector_angle: pi / 2 - incidence_angle,
    projection_law.vector_length: outer_travel_distance,
    projection_law.projection: x
})
projection_incidence_distance = solve(projection_incidence_eq, outer_travel_distance)[0]
outer_travel_distance_polar_eq = Eq(outer_travel_distance, projection_incidence_distance)
projection_refraction_eq = projection_law.law.subs({
    projection_law.vector_angle: pi / 2 - refraction_angle,
    projection_law.vector_length: medium_travel_distance,
    projection_law.projection: d - x
})
projection_refraction_distance = solve(projection_refraction_eq, medium_travel_distance)[0]
medium_travel_distance_polar_eq = Eq(medium_travel_distance, projection_refraction_distance)

# NOTE: derivative of the path equals to zero does not always correspond to minumum travel time. It
#       can also be a maximum time. Here we prove that Snell's law is applicable to cases
#       with extreme travel time.
# NOTE: light travel path is known to have minimum time, which is an extreme time case.

# Substitute l1, l2 in travel_time, differentiate by x and equate the derivative to zero - extreme time case.
travel_time_on_cartesian = travel_time.subs({
    outer_travel_distance_cartesian_eq.lhs: outer_travel_distance_cartesian_eq.rhs,
    medium_travel_distance_cartesian_eq.lhs: medium_travel_distance_cartesian_eq.rhs
})
extreme_time_case = Eq(diff(travel_time_on_cartesian, x), 0)

# Let's get the same expression using the sines of the angles.
extreme_time_case = extreme_time_case.subs({
    outer_travel_distance_cartesian_eq.rhs: outer_travel_distance_cartesian_eq.lhs,
    medium_travel_distance_cartesian_eq.rhs: medium_travel_distance_cartesian_eq.lhs
})
extreme_time_case = extreme_time_case.subs({
    outer_travel_distance_polar_eq.lhs: outer_travel_distance_polar_eq.rhs,
    medium_travel_distance_polar_eq.lhs: medium_travel_distance_polar_eq.rhs
}).simplify()

# Finally, let's use the definition of the refractive index as the ratio of the speed of light in the medium to that in a reference medium (vacuum).
outer_refraction_definition = refractive_index_definition.definition.subs(
    {refractive_index_definition.relative_refractive_index: incidence_refractive_index})
medium_refreaction_definition = refractive_index_definition.definition.subs(
    {refractive_index_definition.relative_refractive_index: resulting_refractive_index})
outer_refraction_velocity = solve(outer_refraction_definition,
    refractive_index_definition.refracted_wave_speed,
    dict=True)[0][refractive_index_definition.refracted_wave_speed]
medium_refraction_velocity = solve(medium_refreaction_definition,
    refractive_index_definition.refracted_wave_speed,
    dict=True)[0][refractive_index_definition.refracted_wave_speed]

extreme_time_case = extreme_time_case.subs({
    outer_velocity: outer_refraction_velocity,
    medium_velocity: medium_refraction_velocity
})
resulting_expression = Eq(
    incidence_refractive_index * sin(incidence_angle),
    solve(extreme_time_case, incidence_refractive_index * sin(incidence_angle),
    dict=True)[0][incidence_refractive_index * sin(incidence_angle)])

# Verify that the resulting_expression corresponds to the law.
assert expr_equals(law.lhs, resulting_expression.lhs)
assert expr_equals(law.rhs, resulting_expression.rhs)


def print_law() -> str:
    return print_expression(law)


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
