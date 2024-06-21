from sympy import Expr, symbols
from symplyphysics import (
    units,
    angle_type,
    validate_input,
    validate_output,
    Vector,
    Quantity,
    QuantityVector,
    scale_vector,
    add_cartesian_vectors,
    cross_cartesian_vectors,
    subtract_cartesian_vectors,
)
from symplyphysics.core.vectors.arithmetics import diff_cartesian_vector

# Description
## Suppose two reference frames, one of which is fixed (S) and the other one is moving (S'). The movement of
## a body stationary in moving frame S' due to the movement of the frame itself is called transfer movement.
## The acceleration related to such movement is called transfer acceleration. It is composed of the acceleration
## of the moving frame relative to the fixed frame, centripetal acceleration and the acceleration due to uneven
## rotation of the moving frame. The transfer acceleration only depends on the movement of frame S' relative to
## stationary frame S, so it would be the acceleration in S of a point stationary in S'.

# Law: a_tr = a_0 + cross(w, cross(w, r)) + cross(d(w)/dt, r)
## a_tr - vector of transfer acceleration of point X relative to fixed frame S
## a_0 - vector of acceleration of moving frame S' relative to fixed frame S
## w - pseudovector of angular velocity related to rotation of moving frame S' about instantaneous axis
## r - vector of position of point X relative to moving frame S'
## d/dt - derivative w.r.t. time


def transfer_acceleration_law(
    moving_frame_acceleration_: Vector,
    angular_velocity_: Vector,
    position_vector_: Vector,
    time_: Expr,
) -> Vector:
    centripetal_acceleration_ = cross_cartesian_vectors(
        angular_velocity_,
        cross_cartesian_vectors(angular_velocity_, position_vector_),
    )
    rotation_acceleration = cross_cartesian_vectors(
        diff_cartesian_vector(angular_velocity_, time_),
        position_vector_,
    )
    return add_cartesian_vectors(
        moving_frame_acceleration_,
        add_cartesian_vectors(
            centripetal_acceleration_,
            rotation_acceleration,
        )
    )


@validate_input(
    moving_frame_acceleration_=units.acceleration,
    angular_velocity_before_=angle_type / units.time,
    angular_velocity_after_=angle_type / units.time,
    position_vector_=units.length,
    time_before_=units.time,
    time_after_=units.time,
    time_=units.time,
)
@validate_output(units.acceleration)
def calculate_transfer_acceleration(
    moving_frame_acceleration_: QuantityVector,
    angular_velocity_before_: QuantityVector,
    angular_velocity_after_: QuantityVector,
    position_vector_: QuantityVector,
    time_before_: Quantity,
    time_after_: Quantity,
    time_: Quantity,
) -> QuantityVector:
    time = symbols("time")

    # (w - w0) / (t - t0) = (w1 - w0) / (t1 - t0)
    angular_velocity_ = add_cartesian_vectors(
        angular_velocity_before_.to_base_vector(),
        scale_vector(
            (time - time_before_) / (time_after_ - time_before_),
            subtract_cartesian_vectors(
                angular_velocity_after_.to_base_vector(),
                angular_velocity_before_.to_base_vector(),
            ),
        ),
    )

    result = transfer_acceleration_law(
        moving_frame_acceleration_.to_base_vector(),
        angular_velocity_,
        position_vector_.to_base_vector(),
        time,
    )

    return QuantityVector.from_base_vector(
        result,
        subs={time: time_},
    )
