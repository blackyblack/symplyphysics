from pytest import approx
from symplyphysics import (
    Quantity,
    dot_vectors,
    units,
    angle_type,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    cross_cartesian_vectors,
)

# Description
## Assuming a body rotating about a fixed axis, the vector of its linear displacement can be expressed
## as the cross product of the vector of angular displacement and the radius of rotation.

# Law: s = cross(theta, r)
## s - vector of linear displacement
## theta - pseudovector of angular displacement, parallel to axis of rotation
## r - radius vector of body, perpendicular to axis of rotation
## cross(a, b) - cross product between vectors a and b

# Condition
## Angular displacement pseudovector and radius vector should be perpendicular to each other


def linear_displacement_law(angular_displacement_: Vector, rotation_radius_: Vector) -> Vector:
    return cross_cartesian_vectors(angular_displacement_, rotation_radius_)


@validate_input(angular_displacement_=angle_type, rotation_radius_=units.length)
@validate_output(units.length)
def calculate_linear_displacement(angular_displacement_: QuantityVector,
    rotation_radius_: QuantityVector) -> QuantityVector:
    angular_displacement_vector = angular_displacement_.to_base_vector()
    rotation_radius_vector = rotation_radius_.to_base_vector()
    dot_vectors_result = Quantity(dot_vectors(angular_displacement_vector, rotation_radius_vector))
    if dot_vectors_result.scale_factor != approx(0.0, rel=1e-3):
        raise ValueError(
            "Angular displacement pseudovector and rotation radius vector should be perpendicular to each other"
        )
    result_vector = linear_displacement_law(angular_displacement_vector, rotation_radius_vector)
    return QuantityVector.from_base_vector(result_vector)
