from symplyphysics import (
    units,
    Symbol,
    Vector,
    Quantity,
    QuantityVector,
    scale_vector,
    validate_input,
    validate_output,
)
from symplyphysics.core.symbols.quantities import list_of_quantities

# Description
## A particle's linear momentum is a vector quantity defines as its velocity vector multiplied by its mass.

# Definition: p = m * v
## p - vector of linear momentum
## m - mass
## v - velocity vector

mass = Symbol("mass", units.mass)


def momentum_definition(velocity_: Vector) -> Vector:
    return scale_vector(mass, velocity_)


def velocity_from_momentum_law(momentum_: Vector) -> Vector:
    return scale_vector(1 / mass, momentum_)


@validate_input(mass_=mass, velocity_=units.velocity)
@validate_output(units.momentum)
def calculate_momentum(mass_: Quantity, velocity_: QuantityVector) -> QuantityVector:
    momentum = momentum_definition(velocity_)
    momentum_components = list_of_quantities(momentum.components, {mass: mass_})
    return QuantityVector(momentum_components, velocity_.coordinate_system)


@validate_input(mass_=mass, momentum_=units.momentum)
@validate_output(units.velocity)
def calculate_velocity_from_momentum(mass_: Quantity, momentum_: QuantityVector) -> QuantityVector:
    velocity = velocity_from_momentum_law(momentum_)
    velocity_components = list_of_quantities(velocity.components, {mass: mass_})
    return QuantityVector(velocity_components, momentum_.coordinate_system)
