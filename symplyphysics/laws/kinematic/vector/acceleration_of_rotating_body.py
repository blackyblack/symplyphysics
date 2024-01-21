from symplyphysics import (
    units,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    add_cartesian_vectors,
)

# Description
## The acceleration of a rotating body is composed of two parts: radial acceleration, which is always
## present in a rotating environment and which points to the center of rotation, and tangential acceleration,
## which is responsible for the change of the magnitude of the velocity vector.

# Law: a = a_r + a_t
## a - total acceleration vector
## a_r - radial acceleration vector (parallel to position vector of the rotating point)
## a_t - tangential acceleration vector (perpendicular to position vector of the rotating
##       point, and tangent to its path)


def acceleration_law(
    radial_acceleration_: Vector,
    tangential_acceleration_: Vector,
) -> Vector:
    return add_cartesian_vectors(radial_acceleration_, tangential_acceleration_)


@validate_input(
    radial_acceleration_=units.acceleration,
    tangential_acceleration_=units.acceleration,
)
@validate_output(units.acceleration)
def calculate_acceleration(
    radial_acceleration_: QuantityVector,
    tangential_acceleration_: QuantityVector,
) -> QuantityVector:
    acceleration_vector = acceleration_law(
        radial_acceleration_,
        tangential_acceleration_,
    )
    return QuantityVector(
        acceleration_vector.components,
        radial_acceleration_.coordinate_system,
    )
