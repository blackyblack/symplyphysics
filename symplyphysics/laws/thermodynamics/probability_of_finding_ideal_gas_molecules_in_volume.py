from sympy import Eq
from symplyphysics import (
    Quantity,
    Symbol,
    dimensionless,
    units,
    validate_input,
    validate_output,
    convert_to_float,
)
from symplyphysics.core.symbols.probability import Probability

# Description
## For ideal gas molecules, we can calculate the probability of N molecules occupying a certain region
## in space, say, a container with gas particles.

# Law: P = (V / V0)**N
## P - probability
## V0 - total volume of container
## V - volume of part of container
## N - particle count

# Conditions
## - Gas is ideal
## - `V <= V0`
## - There are no external fields acting on the particles

probability = Symbol("probability", dimensionless)
total_volume = Symbol("total_volume", units.volume)
partial_volume = Symbol("partial_volume", units.volume)
particle_count = Symbol("particle_count", dimensionless, integer=True)

law = Eq(probability, (partial_volume / total_volume)**particle_count)


@validate_input(
    total_volume_=total_volume,
    partial_volume_=partial_volume,
    particle_count_=particle_count,
)
@validate_output(probability)
def calculate_probability(
    total_volume_: Quantity,
    partial_volume_: Quantity,
    particle_count_: int,
) -> Probability:
    if partial_volume_.scale_factor > total_volume_.scale_factor:
        raise ValueError("Partial volume must not exceed total volume")

    result = law.rhs.subs({
        total_volume: total_volume_,
        partial_volume: partial_volume_,
        particle_count: particle_count_,
    })

    return Probability(convert_to_float(result))
