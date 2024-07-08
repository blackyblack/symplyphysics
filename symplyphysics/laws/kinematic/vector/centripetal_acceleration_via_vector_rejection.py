from symplyphysics import (
    units,
    angle_type,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    vector_magnitude,
    scale_vector,
)
from symplyphysics.core.vectors.arithmetics import reject_cartesian_vector


# Description
## Centripetal acceleration is the acceleration of a body in a rotating coordinate system
## which is directed towards the axis of rotation.

# Law: a_c = -1 * norm(w)**2 * reject(r, w)
## a_c - vector of centripetal acceleration
## w - pseudovector of angular velocity
## r - vector of position of body
## norm(a) - Euclidean norm of vector a
## reject(a, b) - component of vector a perpendicular to vector b, i.e. rejection of a from b


def centripetal_acceleration_law(
    angular_velocity_: Vector,
    position_vector_: Vector,
) -> Vector:
    return scale_vector(
        -1 * vector_magnitude(angular_velocity_)**2,
        reject_cartesian_vector(position_vector_, angular_velocity_),
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
