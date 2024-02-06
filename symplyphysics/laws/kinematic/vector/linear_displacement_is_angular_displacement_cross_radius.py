from pytest import approx
from symplyphysics import (
    units,
    angle_type,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    cross_cartesian_vectors,
)
from symplyphysics.core.vectors.arithmetics import dot_quantity_vectors

# Description
## Assuming a body rotating about a fixed axis, the vector of its linear displacement can be expressed
## as the cross product of the vector of angular displacement and the radius of rotation.

# Law: s = [theta, r]
## s - vector of linear displacement
## theta - pseudovector of angular displacement, parallel to axis of rotation
## r - radius vector of body, perpendicular to axis of rotation
## [a, b] - cross product between vectors a and b

# Condition
## Angular displacement pseudovector and radius vector should be perpendicular to each other


def linear_displacement_law(angular_displacement_: Vector, rotation_radius_: Vector) -> Vector:
    return cross_cartesian_vectors(angular_displacement_, rotation_radius_)


@validate_input(angular_displacement_=angle_type, rotation_radius_=units.length)
@validate_output(units.length)
def calculate_linear_displacement(angular_displacement_: QuantityVector,
    rotation_radius_: QuantityVector) -> QuantityVector:
    if dot_quantity_vectors(angular_displacement_, rotation_radius_).scale_factor != approx(0.0,
        rel=1e-3):
        raise ValueError(
            "Angular displacement psedovector and rotation radius vector should be perpendular to each other."
        )
    result_vector = linear_displacement_law(angular_displacement_, rotation_radius_)
    return QuantityVector(result_vector.components, angular_displacement_.coordinate_system)
