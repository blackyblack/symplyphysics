from symplyphysics import (
    units,
    Vector,
    QuantityVector,
    cross_cartesian_vectors,
    validate_input,
    validate_output,
)

# Description
## The angular momentum of a particle is defined relative to a fixed point, usually an origin,
## as a cross product of its position vector and linear momentum.

# Definition: L = cross(r, p)
## L - vector of angular momentum
## p - vector of linear momentum
## r - position vector of particle relative to a fixed point
## cross(a, b) - cross product between vectors a and b


def angular_momentum_law(position_vector_: Vector, linear_momentum_: Vector) -> Vector:
    return cross_cartesian_vectors(position_vector_, linear_momentum_)


@validate_input(position_vector_=units.length, linear_momentum_=units.momentum)
@validate_output(units.length * units.momentum)
def calculate_angular_momentum(
    position_vector_: QuantityVector, linear_momentum_: QuantityVector
) -> QuantityVector:
    result = angular_momentum_law(position_vector_, linear_momentum_)
    return QuantityVector(result.components, position_vector_.coordinate_system)
