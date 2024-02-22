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

# Description
## A particle's linear momentum is a vector quantity defined as its velocity vector multiplied by its mass.

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
    result_vector = momentum_definition(velocity_.to_base_vector())
    return QuantityVector.from_base_vector(
        result_vector, subs={mass: mass_}
    )


@validate_input(mass_=mass, momentum_=units.momentum)
@validate_output(units.velocity)
def calculate_velocity(mass_: Quantity, momentum_: QuantityVector) -> QuantityVector:
    result_vector = velocity_from_momentum_law(momentum_.to_base_vector())
    return QuantityVector.from_base_vector(
        result_vector, subs={mass: mass_}
    )
