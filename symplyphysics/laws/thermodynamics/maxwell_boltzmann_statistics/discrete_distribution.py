from sympy import Eq, exp, solve
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
)

# Description
## Maxwell-Boltzmann distribution can be written as a discrete distribution of a single particle's
## discrete energy spectrum. Maxwell-Boltzmann statistics gives the average number of particles
## found in a given single-particle microstate.

# Law: N_i / N = exp(-E_i / (k * T)) / Z
## N_i - expected number of particles in the single-particle microstate i
## N - total number of particles in the system
## E_i - energy of microstate i
## k - Boltzmann constant
## T - equilibrium temperature of the system
## Z - single-particle partition function (normalizing factor)

# Conditions
## - Particles do not interact and are classical
## - The system is in thermal equilibrium.

particle_count_in_microstate = Symbol("particle_count_in_microstate", dimensionless)
total_particle_count = Symbol("total_particle_count", dimensionless)
energy_of_microstate = Symbol("energy_of_microstate", units.energy)
equilibrium_temperature = clone_symbol(symbols.thermodynamics.temperature, "equilibrium_temperature")
single_particle_partition_function = Symbol("single_particle_partition_function", dimensionless)

law = Eq(
    particle_count_in_microstate / total_particle_count,
    exp(-1 * energy_of_microstate / (units.boltzmann_constant * equilibrium_temperature))
    / single_particle_partition_function
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    total_particle_count_=total_particle_count,
    energy_of_microstate_=energy_of_microstate,
    equilibrium_temperature_=equilibrium_temperature,
    single_particle_partition_function_=single_particle_partition_function,
)
@validate_output(particle_count_in_microstate)
def calculate_particle_count_in_microstate(
    total_particle_count_: int,
    energy_of_microstate_: Quantity,
    equilibrium_temperature_: Quantity,
    single_particle_partition_function_: float,
) -> Quantity:
    expr = solve(law, particle_count_in_microstate)[0]

    result = expr.subs({
        total_particle_count: total_particle_count_,
        energy_of_microstate: energy_of_microstate_,
        equilibrium_temperature: equilibrium_temperature_,
        single_particle_partition_function: single_particle_partition_function_,
    })
    return Quantity(result)
