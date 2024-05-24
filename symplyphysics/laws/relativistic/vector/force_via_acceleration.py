from sympy import sqrt
from sympy.physics.units import speed_of_light
from symplyphysics import (
    units,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
    Vector,
    add_cartesian_vectors,
    dot_vectors,
    scale_vector,
    QuantityVector,
)
from symplyphysics.core.vectors.arithmetics import decompose_into_projections

# Description
## In special relativity, the Newton's second law does not hold in the form `F = m * a`. There still exists a relation
## between force `F` and acceleration `a`, albeit more complex.

# Law: F = gamma**3 * m0 * a_parallel + gamma * m0 * a_orthogonal
## F - force vector
## a - acceleration vector
## v - velocity vector
## a_parallel = dot(a, v) / dot(v, v) * v - acceleration component parallel to velocity
## a_orthogonal = a - a_parallel - acceleration component orthogonal to velocity
## m0 - rest mass
## gamma = 1 / sqrt(1 - dot(v, v) / c**2) - Lorentz factor
## dot(a, b) - dot product between vectors `a` and `b`

# Conditions
## - This law applies to special relativity.

rest_mass = Symbol("rest_mass", units.mass)


def force_law(acceleration_: Vector, velocity_: Vector) -> Vector:
    acceleration_parallel_, acceleration_orthogonal_ = decompose_into_projections(
        acceleration_, velocity_
    )

    lorentz_factor_ = 1 / sqrt(1 - dot_vectors(velocity_, velocity_) / speed_of_light**2)

    force_parallel_ = scale_vector(
        lorentz_factor_**3 * rest_mass,
        acceleration_parallel_,
    )
    force_orthogonal_ = scale_vector(
        lorentz_factor_ * rest_mass,
        acceleration_orthogonal_,
    )

    return add_cartesian_vectors(
        force_parallel_,
        force_orthogonal_,
    )


@validate_input(
    rest_mass_=rest_mass,
    acceleration_=units.acceleration,
    velocity_=units.velocity,
)
@validate_output(units.force)
def calculate_force(
    rest_mass_: Quantity,
    acceleration_: QuantityVector,
    velocity_: QuantityVector,
) -> QuantityVector:
    result = force_law(
        acceleration_.to_base_vector(),
        velocity_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(
        result,
        subs={rest_mass: rest_mass_},
    )
