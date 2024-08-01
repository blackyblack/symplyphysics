from sympy import Expr, abc
from symplyphysics import (
    units,
    angle_type,
    validate_input,
    validate_output,
    Vector,
    Quantity,
    QuantityVector,
    cross_cartesian_vectors,
    scale_vector,
)
from symplyphysics.core.vectors.arithmetics import diff_cartesian_vector

# Description
## Suppose two reference frames, one of which is fixed (S) and the other is moving (S'). When S' rotates
## around S in a non-uniform way, the acceleration of some body B in S has a component corresponding to
## that non-uniform rotation of S'. It is part of the transfer acceleration of body B in S.

# Law: a_rot = cross(dw/dt, r)
## a_rot - vector of body B's acceleration due to non-uniform rotation of S'
## w - pseudovector of angular velocity of rotation of S'
## d/dt - time derivative
## r - position vector of body B within S'
## cross(a, b) - vector of cross product between vectors a and b


def non_uniform_rotation_acceleration_law(
    angular_velocity_: Vector,
    time_: Expr,
    position_vector_: Vector,
) -> Vector:
    return cross_cartesian_vectors(
        diff_cartesian_vector(angular_velocity_, time_),
        position_vector_,
    )


@validate_input(
    angular_velocity_change_=angle_type / units.time,
    time_change_=units.time,
    position_vector=units.length,
)
@validate_output(units.acceleration)
def calculate_non_uniform_rotation_acceleration(
    angular_velocity_change_: QuantityVector,
    time_change_: Quantity,
    position_vector_: QuantityVector,
) -> QuantityVector:
    time_ = abc.t

    angular_velocity_ = scale_vector(
        time_ / time_change_,
        angular_velocity_change_.to_base_vector(),
    )

    result_vector = non_uniform_rotation_acceleration_law(
        angular_velocity_,
        time_,
        position_vector_.to_base_vector(),
    )

    return QuantityVector.from_base_vector(result_vector)
