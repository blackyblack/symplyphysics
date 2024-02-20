from pytest import approx
from symplyphysics import (
    Quantity,
    dot_vectors,
    units,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    angle_type,
    cross_cartesian_vectors,
)

# Description
## Assuming the body is rotating about a fixed axis, the vector of its linear velocity is the
## cross product of the vector of angular rotation and the radius vector of the body. The radius
## vector of the body is perpendicular to the axis of rotation.

# Law: v = cross(w, r)
## v - vector of linear velocity
## w - pseudovector of angular velocity, parallel to axis of rotation
## r - radius vector of body, perpendicular to axis of rotation
## cross(a, b) - cross product between vectors a and b

# Conditions:
## - Angular velocity pseudovector and radius vector should be perpendicular to each other


def linear_velocity_law(angular_velocity: Vector, rotation_radius: Vector) -> Vector:
    return cross_cartesian_vectors(angular_velocity, rotation_radius)


@validate_input(angular_velocity_=angle_type / units.time, rotation_radius_=units.length)
@validate_output(units.velocity)
def calculate_linear_velocity(angular_velocity_: QuantityVector,
    rotation_radius_: QuantityVector) -> QuantityVector:
    angular_velocity_vector = angular_velocity_.to_base_vector()
    rotation_radius_vector = rotation_radius_.to_base_vector()
    dot_vectors_result = Quantity(dot_vectors(angular_velocity_vector, rotation_radius_vector))
    if dot_vectors_result.scale_factor != approx(0.0, rel=1e-3):
        raise ValueError(
            "Angular velocity pseudovector and rotation radius vector should be perpendicular to each other"
        )
    result_vector = linear_velocity_law(angular_velocity_vector, rotation_radius_vector)
    return QuantityVector.from_base_vector(result_vector)
