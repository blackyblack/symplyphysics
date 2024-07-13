from sympy import Expr, Rational
from symplyphysics import (
    angle_type,
    units,
    Quantity,
    Vector,
    QuantityVector,
    validate_input,
    validate_output,
    scale_vector,
    add_cartesian_vectors,
    cross_cartesian_vectors,
)

# Description
## Suppose a reference frame S' is fixed to a moving object A (e.g., Earth) and some body B moving freely
## (i.e. the sum of external non-gravitational forces acting on it is zero). In the case of an inertial
## frame of reference, the displacement of body B from the starting position would follow the usual rule
## `s = v0 * t + (g / 2) * t**2`. But in the case of non-inertial frames, we have to take the Coriolis and
## the centrifugal force into account as well, which results into the following series:

# Conditions
## - The sum `F` of all other, non-gravitational forces acting on body B is 0.
## - `g` is independent of coordinates.

# Notes
## - Note that the series is truncated at the fifth term. More terms can be obtained by plugging the result
##   into the equation of motion `a = g + cross(v, w)` and integrating it over time.

# Law: s = v0 * t
#        + t**2 * (g / 2 + cross(v0, w))
#        + t**3 / 3 * (cross(g, w) + 2 * cross(cross(v0, w), w))
#        + t**4 / 6 * cross(cross(g, w), w)
#        + O(t**5)
## t - time
## s - vector of displacement of body B in S'
## v0 - vector of initial velocity of body B
## w - pseudovector of angular velocity of rotation of moving frame S'
## g - vector of acceleration due to gravity (or free fall acceleration) of body B
## cross(a, b) - cross product between vectors `a` and `b`
## O - [big O](https://en.wikipedia.org/wiki/Big_O_notation) as per the formal mathematical definition
##     in the limit of t approaching infinity


def displacement_law(
    time_: Expr,
    initial_velocity_: Vector,
    angular_velocity_: Vector,
    acceleration_due_to_gravity_: Vector,
) -> Vector:
    initial_cross_angular = cross_cartesian_vectors(
        initial_velocity_,
        angular_velocity_,
    )

    acceleration_cross_angular = cross_cartesian_vectors(
        acceleration_due_to_gravity_,
        angular_velocity_,
    )

    return add_cartesian_vectors(
        scale_vector(time_, initial_velocity_),
        scale_vector(
        time_**2,
        add_cartesian_vectors(
        scale_vector(Rational(1, 2), acceleration_due_to_gravity_),
        initial_cross_angular,
        ),
        ),
        scale_vector(
        time_**3 / 3,
        add_cartesian_vectors(
        acceleration_cross_angular,
        scale_vector(2, cross_cartesian_vectors(initial_cross_angular, angular_velocity_)),
        ),
        ),
        scale_vector(time_**4 / 6,
        cross_cartesian_vectors(acceleration_cross_angular, angular_velocity_)))


@validate_input(
    time_=units.time,
    initial_velocity_=units.velocity,
    angular_velocity_=angle_type / units.time,
    acceleration_due_to_gravity_=units.acceleration,
)
@validate_output(units.length)
def calculate_displacement(
    time_: Quantity,
    initial_velocity_: QuantityVector,
    angular_velocity_: QuantityVector,
    acceleration_due_to_gravity_: QuantityVector,
) -> QuantityVector:
    result_vector = displacement_law(
        time_,
        initial_velocity_.to_base_vector(),
        angular_velocity_.to_base_vector(),
        acceleration_due_to_gravity_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(result_vector)
