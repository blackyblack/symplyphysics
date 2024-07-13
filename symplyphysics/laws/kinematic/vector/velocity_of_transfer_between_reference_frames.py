from symplyphysics import (
    units,
    angle_type,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    add_cartesian_vectors,
    scale_vector,
    cross_cartesian_vectors,
)

# Description
## Suppose two reference frames, one of which is fixed (S) and the other one is moving (S'). The movement of
## a body stationary in moving frame S' due to the movement of the frame itself is called transfer movement.
## The velocity related to such movement is called transfer velocity. For any material point X, its transfer
## velocity relative to fixed frame S is the sum of the velocity of frame S' relative to frame S
## and the cross product of the angular velocity of moving frame's rotation and the position
## vector of X in moving frame S'.

# Law: v_tr = v_0 + cross(w, r)
## v_tr - vector of transfer velocity of point X relative to fixed frame S.
## v_0 - vector of moving frame S' relative to fixed frame S
## w - pseudovector of angular velocity related to rotation of moving frame S' about instantaneous axis
## r - vector of position of point X relative to moving frame S'


def transfer_velocity_law(
    moving_frame_velocity_: Vector,
    angular_velocity_: Vector,
    position_vector_: Vector,
) -> Vector:
    return add_cartesian_vectors(
        moving_frame_velocity_,
        cross_cartesian_vectors(angular_velocity_, position_vector_),
    )


def moving_frame_velocity_law(
    transfer_velocity_: Vector,
    angular_velocity_: Vector,
    position_vector_: Vector,
) -> Vector:
    return add_cartesian_vectors(
        transfer_velocity_,
        scale_vector(-1, cross_cartesian_vectors(angular_velocity_, position_vector_)))


@validate_input(
    moving_frame_velocity_=units.velocity,
    angular_velocity_=angle_type / units.time,
    position_vector_=units.length,
)
@validate_output(units.velocity)
def calculate_transfer_velocity(
    moving_frame_velocity_: QuantityVector,
    angular_velocity_: QuantityVector,
    position_vector_: QuantityVector,
) -> QuantityVector:
    result = transfer_velocity_law(
        moving_frame_velocity_.to_base_vector(),
        angular_velocity_.to_base_vector(),
        position_vector_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(result)
