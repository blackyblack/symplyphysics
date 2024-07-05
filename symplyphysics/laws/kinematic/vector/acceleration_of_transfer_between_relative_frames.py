from symplyphysics import (
    units,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    add_cartesian_vectors,
    subtract_cartesian_vectors,
)

# Description
## Suppose two reference frames, one of which is fixed (S) and the other one is moving (S'). The movement of
## a body stationary in moving frame S' due to the movement of the frame itself is called transfer movement.
## The acceleration related to such movement is called transfer acceleration. It is composed of the acceleration
## of the moving frame relative to the fixed frame, centripetal acceleration and the acceleration due to uneven
## rotation of the moving frame. The transfer acceleration only depends on the movement of frame S' relative to
## stationary frame S, so it would be the acceleration in S of a point stationary in S'.

# Law: a_tr = a_0 + a_centripetal + a_non_uniform_rotation
## a_tr - vector of transfer acceleration of body B relative to fixed frame S
## a_0 - vector of acceleration of moving frame S' relative to fixed frame S
## a_centripetal - vector of centripetal acceleration of body B in moving frame S'
## a_non_uniform_rotation - vector of acceleration caused by non-uniform rotation of frame S'


def transfer_acceleration_law(
    moving_frame_acceleration_: Vector,
    centripetal_acceleration_: Vector,
    rotation_acceleration_: Vector,
) -> Vector:
    return add_cartesian_vectors(
        moving_frame_acceleration_,
        add_cartesian_vectors(
            centripetal_acceleration_,
            rotation_acceleration_,
        )
    )


# a_0 = a_tr - (a_centripetal + a_non_uniform_rotation)
def moving_frame_acceleration_law(
    transfer_acceleration_: Vector,
    centripetal_acceleration_: Vector,
    rotation_acceleration_: Vector,
) -> Vector:
    return subtract_cartesian_vectors(
        transfer_acceleration_,
        add_cartesian_vectors(
            centripetal_acceleration_,
            rotation_acceleration_,
        )
    )


# a_centripetal = a_tr - (a_0 + a_non_uniform_rotation)
def centripetal_acceleration_law(
    transfer_acceleration_: Vector,
    moving_frame_acceleration_: Vector,
    rotation_acceleration_: Vector,
) -> Vector:
    return subtract_cartesian_vectors(
        transfer_acceleration_,
        add_cartesian_vectors(
            moving_frame_acceleration_,
            rotation_acceleration_,
        )
    )


# a_non_uniform_rotation = a_tr - (a_0 + a_centripetal)
def rotation_acceleration_law(
    transfer_acceleration_: Vector,
    moving_frame_acceleration_: Vector,
    centripetal_acceleration_: Vector,
) -> Vector:
    return subtract_cartesian_vectors(
        transfer_acceleration_,
        add_cartesian_vectors(
            moving_frame_acceleration_,
            centripetal_acceleration_,
        )
    )


@validate_input(
    moving_frame_acceleration_=units.acceleration,
    centripetal_acceleration_=units.acceleration,
    rotation_acceleration_=units.acceleration,
)
@validate_output(units.acceleration)
def calculate_transfer_acceleration(
    moving_frame_acceleration_: QuantityVector,
    centripetal_acceleration_: QuantityVector,
    rotation_acceleration_: QuantityVector,
) -> QuantityVector:
    vector = transfer_acceleration_law(
        moving_frame_acceleration_.to_base_vector(),
        centripetal_acceleration_.to_base_vector(),
        rotation_acceleration_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(vector)
