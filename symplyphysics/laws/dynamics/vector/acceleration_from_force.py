from symplyphysics import (units, Quantity, Symbol, QuantityVector, Vector, scale_vector,
    validate_input, validate_output, list_of_quantities)

# Description
## Newton's second law in vector form: a = 1/m * F
## Where:
## F - force vector,
## m - mass,
## a - acceleration vector,
## * - scalar multiplication (scale vector).

mass = Symbol("mass", units.mass)


def acceleration_law(force_: Vector) -> Vector:
    return scale_vector(1 / mass, force_)


def force_law(acceleration_: Vector) -> Vector:
    return scale_vector(mass, acceleration_)


@validate_input(mass_=mass, acceleration_=units.acceleration)
@validate_output(units.force)
def calculate_force(mass_: Quantity, acceleration_: QuantityVector) -> QuantityVector:
    result_force = force_law(acceleration_)
    force_components = list_of_quantities(result_force.components, {mass: mass_})
    return QuantityVector(force_components, acceleration_.coordinate_system)


@validate_input(mass_=mass, force_=units.force)
@validate_output(units.acceleration)
def calculate_acceleration(mass_: Quantity, force_: QuantityVector) -> QuantityVector:
    result_acceleration = acceleration_law(force_)
    acceleration_components = list_of_quantities(result_acceleration.components, {mass: mass_})
    return QuantityVector(acceleration_components, force_.coordinate_system)
