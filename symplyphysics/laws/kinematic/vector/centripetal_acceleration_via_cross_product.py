from symplyphysics import (
    units,
    angle_type,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    cross_cartesian_vectors,
)


# Description
## Centripetal acceleration is the acceleration of a body in a rotating coordinate system
## which is directed towards the axis of rotation.

# Law: a_c = cross(w, cross(w, r))
## a_c - vector of centripetal acceleration
## w - pseudovector of angular velocity
## r - vector of position of body


def centripetal_acceleration_law(
    angular_velocity_: Vector,
    position_vector_: Vector,
) -> Vector:
    return cross_cartesian_vectors(
        angular_velocity_,
        cross_cartesian_vectors(angular_velocity_, position_vector_),
    )


@validate_input(
    angular_velocity_=angle_type / units.time,
    position_vector_=units.length,
)
@validate_output(units.acceleration)
def calculate_centripetal_acceleration(
    angular_velocity_: QuantityVector,
    position_vector_: QuantityVector,
) -> QuantityVector:
    vector = centripetal_acceleration_law(
        angular_velocity_.to_base_vector(),
        position_vector_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(vector)
