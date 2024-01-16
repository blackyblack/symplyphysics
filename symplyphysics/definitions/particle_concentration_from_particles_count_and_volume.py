from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression,
                           validate_input,
                           validate_output, dimensionless)

# Description
## The particle concentration is a physical quantity equal to the ratio of the number of particles
## to the volume in which they are located

# Definition: n = N / V
# Where:
## N is the count of particles
## V is volume
## n is the particle concentration

count_of_particles = Symbol("count_of_particles", dimensionless)
volume = Symbol("volume", units.volume)
concentration = Symbol("concentration", dimensionless / units.volume)

definition = Eq(concentration, count_of_particles / volume)

definition_units_SI = units.meter**(-3)


def print_law() -> str:
    return print_expression(definition)


@validate_input(count_of_particles_=count_of_particles, volume_=volume)
@validate_output(concentration)
def calculate_concentration(count_of_particles_: Quantity, volume_: Quantity) -> Quantity:
    solved = solve(definition, concentration, dict=True)[0][concentration]
    result_expr = solved.subs({count_of_particles: count_of_particles_, volume: volume_})
    return Quantity(result_expr)
