from pytest import approx
from symplyphysics import (
    units,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    angle_type,
    cross_cartesian_vectors,
)
from symplyphysics.core.vectors.arithmetics import dot_quantity_vectors

# Description
## Assuming the body is rotating about a fixed axis, the vector of its linear velocity is the
## cross product of the vector of angular rotation and the radius vector of the body. The radius
## vector of the body is perpendicular to the axis of rotation.

# Law: v = [w, r]
## v - vector of linear velocity
## w - pseudovector of angular velocity, parallel to axis of rotation
## r - radius vector of body, perpendicular to axis of rotation
## [a, b] - cross product between vectors a and b

# Conditions:
## - Angular velocity pseudovector and radius vector should be perpendicular to each other


def linear_velocity_law(angular_velocity: Vector, rotation_radius: Vector) -> Vector:
    return cross_cartesian_vectors(angular_velocity, rotation_radius)


@validate_input(angular_velocity_=angle_type / units.time, rotation_radius_=units.length)
@validate_output(units.velocity)
def calculate_linear_velocity(angular_velocity_: QuantityVector,
    rotation_radius_: QuantityVector) -> QuantityVector:
    if dot_quantity_vectors(angular_velocity_, rotation_radius_).scale_factor != approx(0.0,
        rel=1e-3):
        raise ValueError(
            "Angular velocity pseudovector and rotation radius vector should be perpendicular to each other."
        )
    result_vector = linear_velocity_law(angular_velocity_, rotation_radius_)
    return QuantityVector(result_vector.components, angular_velocity_.coordinate_system)
