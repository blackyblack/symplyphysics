from pytest import approx
from symplyphysics import (
    CoordinateSystem,
    Quantity,
    dot_vectors,
    units,
    validate_input,
    validate_output,
    Vector,
    QuantityVector,
    add_cartesian_vectors,
    scale_vector,
)

# Description
## The acceleration of a rotating body is composed of two parts: radial acceleration, which is always
## present in a rotating environment and which points to the center of rotation, and tangential acceleration,
## which is responsible for the change of the magnitude of the velocity vector.

# Law: a = a_r + a_t
## a - total acceleration vector
## a_r - radial acceleration vector (parallel to position vector of the rotating point and
##       pointing in the direction opposite to it)
## a_t - tangential acceleration vector (perpendicular to position vector of the rotating
##       point, and tangent to its path)


def acceleration_law(
    radial_acceleration_: Vector,
    tangential_acceleration_: Vector,
) -> Vector:
    return add_cartesian_vectors(radial_acceleration_, tangential_acceleration_)


def radial_acceleration_law(
    total_acceleration_: Vector,
    tangential_acceleration_: Vector,
) -> Vector:
    opposite_tangential_vector = scale_vector(-1, tangential_acceleration_)
    return add_cartesian_vectors(total_acceleration_, opposite_tangential_vector)


def tangential_acceleration_law(
    total_acceleration_: Vector,
    radial_acceleration_: Vector,
) -> Vector:
    opposite_radial_vector = scale_vector(-1, radial_acceleration_)
    return add_cartesian_vectors(total_acceleration_, opposite_radial_vector)


@validate_input(
    radial_acceleration_=units.acceleration,
    tangential_acceleration_=units.acceleration,
)
@validate_output(units.acceleration)
def calculate_acceleration(
    radial_acceleration_: QuantityVector,
    tangential_acceleration_: QuantityVector,
) -> QuantityVector:
    if radial_acceleration_.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        raise ValueError("Radial acceleration vector should be in cartesian coordinate system")
    if tangential_acceleration_.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        raise ValueError("Tangential acceleration vector should be in cartesian coordinate system")
    radial_acceleration_vector = radial_acceleration_.to_base_vector()
    tangential_acceleration_vector = tangential_acceleration_.to_base_vector()
    dot_vectors_result = Quantity(
        dot_vectors(radial_acceleration_vector, tangential_acceleration_vector))
    if dot_vectors_result.scale_factor != approx(0.0, rel=1e-3):
        raise ValueError(
            "Radial and tangential acceleration vectors should be perpendicular to each other")
    acceleration_vector = acceleration_law(
        radial_acceleration_vector,
        tangential_acceleration_vector,
    )
    return QuantityVector.from_base_vector(acceleration_vector)


@validate_input(
    total_acceleration_=units.acceleration,
    tangential_acceleration_=units.acceleration,
)
@validate_output(units.acceleration)
def calculate_radial_acceleration(
    total_acceleration_: QuantityVector,
    tangential_acceleration_: QuantityVector,
) -> QuantityVector:
    if total_acceleration_.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        raise ValueError("Total acceleration vector should be in cartesian coordinate system")
    if tangential_acceleration_.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        raise ValueError("Tangential acceleration vector should be in cartesian coordinate system")
    total_acceleration_vector = total_acceleration_.to_base_vector()
    tangential_acceleration_vector = tangential_acceleration_.to_base_vector()
    radial_acceleration_vector = radial_acceleration_law(
        total_acceleration_vector,
        tangential_acceleration_vector,
    )
    return QuantityVector.from_base_vector(radial_acceleration_vector)


@validate_input(
    total_acceleration_=units.acceleration,
    radial_acceleration_=units.acceleration,
)
@validate_output(units.acceleration)
def calculate_tangential_acceleration(
    total_acceleration_: QuantityVector,
    radial_acceleration_: QuantityVector,
) -> QuantityVector:
    if total_acceleration_.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        raise ValueError("Total acceleration vector should be in cartesian coordinate system")
    if radial_acceleration_.coordinate_system.coord_system_type != CoordinateSystem.System.CARTESIAN:
        raise ValueError("Radial acceleration vector should be in cartesian coordinate system")
    total_acceleration_vector = total_acceleration_.to_base_vector()
    radial_acceleration_vector = radial_acceleration_.to_base_vector()
    tangential_acceleration_vector = tangential_acceleration_law(
        total_acceleration_vector,
        radial_acceleration_vector,
    )
    return QuantityVector.from_base_vector(tangential_acceleration_vector)
