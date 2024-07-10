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
## TODO

# Conditions
## - F = 0
## - g is constant

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
        scale_vector(
            time_**4 / 6,
            cross_cartesian_vectors(acceleration_cross_angular, angular_velocity_)
        )
    )


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
