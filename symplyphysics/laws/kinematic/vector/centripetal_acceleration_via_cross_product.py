from symplyphysics import (
    units,
    angle_type,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    cross_cartesian_vectors,
    vector_magnitude,
    scale_vector,
)
from symplyphysics.core.vectors.arithmetics import reject_cartesian_vector


# Description
## Centripetal acceleration is the acceleration of a body in a rotating coordinate system
## which is directed towards the axis of rotation.

# Law: a_c = cross(w, cross(w, r)) = -1 * norm(w)**2 * r_perp_w
## a_c - vector of centripetal acceleration
## w - pseudovector of angular velocity
## r - vector of position of body
## cross(a, b) - cross product between vectors a and b
## norm(a) - Euclidean norm of vector a
## a_perp_b - component of vector a perpendicular to vector b, i.e. rejection of a from b


def cross_product_law(
    angular_velocity_: Vector,
    position_vector_: Vector,
) -> Vector:
    return cross_cartesian_vectors(
        angular_velocity_,
        cross_cartesian_vectors(angular_velocity_, position_vector_),
    )


def rejection_law(
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
    vector = cross_product_law(
        angular_velocity_.to_base_vector(),
        position_vector_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(vector)
