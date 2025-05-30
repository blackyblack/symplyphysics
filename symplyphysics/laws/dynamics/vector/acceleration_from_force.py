from symplyphysics import (
    units,
    Quantity,
    QuantityVector,
    Vector,
    scale_vector,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

# Description
## Newton's second law in vector form: a = 1/m * F
## Where:
## F - force vector,
## m - mass,
## a - acceleration vector,
## * - scalar multiplication (scale vector).

mass = clone_as_symbol(symbols.mass, positive=True)


def acceleration_law(force_: Vector) -> Vector:
    return scale_vector(1 / mass, force_)


def force_law(acceleration_: Vector) -> Vector:
    return scale_vector(mass, acceleration_)


# TODO: Derive this law from law of force and momentum
# Condition: mass is constant


@validate_input(mass_=mass, acceleration_=units.acceleration)
@validate_output(units.force)
def calculate_force(mass_: Quantity, acceleration_: QuantityVector) -> QuantityVector:
    result_vector = force_law(acceleration_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector, subs={mass: mass_})


@validate_input(mass_=mass, force_=units.force)
@validate_output(units.acceleration)
def calculate_acceleration(mass_: Quantity, force_: QuantityVector) -> QuantityVector:
    result_vector = acceleration_law(force_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector, subs={mass: mass_})
