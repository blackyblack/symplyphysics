from sympy import sqrt
from sympy.physics.units import speed_of_light
from symplyphysics import (
    units,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
    Vector,
    scale_vector,
    QuantityVector,
    add_cartesian_vectors,
    dot_vectors,
)

# Description
# TODO

# Law: a = 1 / (m0 * gamma) * (F - dot(v, F) * v / c**2)
## a - vector of acceleration
## m0 - rest mass
## v - vector of velocity
## c - speed of light
## F - vector of force
## gamma = 1 / sqrt(1 - dot(v, v) / c**2) - Lorentz factor
## dot(a, b) - dot product between vectors `a` and `b`

rest_mass = Symbol("rest_mass", units.mass)


def acceleration_law(force_: Vector, velocity_: Vector) -> Vector:
    force_parallel_to_velocity_ = scale_vector(
        -1 * dot_vectors(velocity_, force_) / speed_of_light**2,
        velocity_,
    )

    resulting_force = add_cartesian_vectors(force_, force_parallel_to_velocity_)

    mass_factor_ = sqrt(1 - dot_vectors(velocity_, velocity_) / speed_of_light**2) / rest_mass

    return scale_vector(mass_factor_, resulting_force)


@validate_input(
    rest_mass_=rest_mass,
    force_=units.force,
    velocity_=units.velocity,
)
@validate_output(units.acceleration)
def calculate_acceleration(
    rest_mass_: Quantity,
    force_: QuantityVector,
    velocity_: QuantityVector,
) -> QuantityVector:
    result = acceleration_law(force_.to_base_vector(), velocity_.to_base_vector())
    return QuantityVector.from_base_vector(
        result,
        subs={rest_mass: rest_mass_},
    )
