from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    Vector,
    add_cartesian_vectors,
    vector_magnitude,
    scale_vector,
    QuantityVector,
    symbols,
)
from symplyphysics.definitions import lorentz_factor as lorentz_factor_def

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

# Links: Wikipedia, end of paragraph <https://en.wikipedia.org/wiki/Relativistic_angular_momentum#Dynamic_mass_moment>

rest_mass = symbols.rest_mass
time = symbols.time


def mass_moment_law(position_: Vector, velocity_: Vector) -> Vector:
    summed_vector_ = add_cartesian_vectors(
        position_,
        scale_vector(-1 * time, velocity_),
    )

    lorentz_factor_ = lorentz_factor_def.definition.rhs.subs({
        lorentz_factor_def.speed: vector_magnitude(velocity_),
    })

    return scale_vector(rest_mass * lorentz_factor_**2, summed_vector_)


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
        subs={
        rest_mass: rest_mass_,
        time: time_
        },
    )
