from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    dimensionless,
    convert_to_float,
    validate_input,
    validate_output,
)

# Description
## The Avogadro constant is the proportionality factor that relates the number of constituent particles
## (usually molecules, atoms or ions) in a sample with the amount of substance in that sample.

## Law: Na = n / mole_count
## Where:
## n is the number of particles
## mole_count is the number of moles of the substance
## Na is Avogadro's number

particles_count = Symbol("particles_count", dimensionless)
mole_count = Symbol("mole_count", units.amount_of_substance)

law = Eq(units.avogadro, particles_count / mole_count)


@validate_input(mole_count_=mole_count)
@validate_output(particles_count)
def calculate_particles_count(mole_count_: Quantity) -> int:
    solved = solve(law, particles_count, dict=True)[0][particles_count]
    result_expr = solved.subs(mole_count, mole_count_)
    result = Quantity(result_expr)
    return int(convert_to_float(result))
