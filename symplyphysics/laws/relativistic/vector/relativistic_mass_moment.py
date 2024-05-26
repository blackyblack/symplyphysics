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

# Description
## Mass moment is an additive physical quantity useful for deriving the Lorentz transformation of angular momentum.
## For isolated systems, it is conserved in time, but unlike angular momentum, the vector of mass moment is a polar
## ("ordinary") vector and therefore invariant under inversion.

# Law: N = m0 * gamma**2 * (x - v * t)
## N - vector of mass moment
## m0 - rest mass
## x - vector of position
## v - vector of velocity
## gamma = 1 / sqrt(1 - dot(v, v)**2 / c**2) - Lorentz factor, where `dot` is dot product and `c` is speed of light
## t - time

rest_mass = Symbol("rest_mass", units.mass)
time = Symbol("time", units.time)


def mass_moment_law(position_: Vector, velocity_: Vector) -> Vector:
    summed_vector = add_cartesian_vectors(
        position_,
        scale_vector(-1 * time, velocity_),
    )

    factor = rest_mass / (1 - dot_vectors(velocity_, velocity_) / speed_of_light**2)

    return scale_vector(factor, summed_vector)


@validate_input(
    rest_mass_=rest_mass,
    position_=units.length,
    velocity_=units.velocity,
    time_=time,
)
@validate_output(units.mass * units.length)
def calculate_mass_moment(
    rest_mass_: Quantity,
    position_: QuantityVector,
    velocity_: QuantityVector,
    time_: Quantity,
) -> QuantityVector:
    result = mass_moment_law(
        position_.to_base_vector(),
        velocity_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(
        result,
        subs={rest_mass: rest_mass_, time: time_},
    )
