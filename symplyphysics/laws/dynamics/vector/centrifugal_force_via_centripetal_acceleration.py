from symplyphysics import (
    units,
    Quantity,
    QuantityVector,
    Vector,
    validate_input,
    validate_output,
    symbols,
    scale_vector,
    vector_magnitude,
)

# Description
## The centrifugal force is a fictitious force that only exists in a non-inertial coordinate system.
## It disappears entirely when transitioning into an inertial frame of reference. It appears to act
## on all objects when viewed in a rotating frame of reference and is directed away from the axis
## of rotation.

# Law: F_cf = -m * a_cp
## F_cf - vector of centrifugal force
## m - mass of body B
## a_cp - vector of centripetal acceleration

mass = symbols.basic.mass


def centrifugal_force_law(centripetal_acceleration_: Vector) -> Vector:
    return scale_vector(-1 * mass, centripetal_acceleration_)


def centripetal_acceleration_law(centrifugal_force_: Vector) -> Vector:
    return scale_vector(-1 / mass, centrifugal_force_)


def mass_law(
    centrifugal_force_: Vector,
    centripetal_acceleration_: Vector,
) -> Vector:
    return vector_magnitude(centrifugal_force_) / vector_magnitude(centripetal_acceleration_)


@validate_input(
    mass_=mass,
    centripetal_acceleration_=units.acceleration,
)
@validate_output(units.force)
def calculate_centrifugal_force(
    mass_: Quantity,
    centripetal_acceleration_: QuantityVector,
) -> QuantityVector:
    result_vector = centrifugal_force_law(centripetal_acceleration_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector, subs={mass: mass_})
