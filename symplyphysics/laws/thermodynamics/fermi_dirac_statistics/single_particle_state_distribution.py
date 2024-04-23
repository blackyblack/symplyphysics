from sympy import Eq, exp
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols,
    convert_to_float,
)

# Description
## Fermi-Dirac statistics is a type of quantum statistics that applies to the physics of a system
## consisting of many non-interacting, identical particles that obey the Pauli exclusion principle.
## For a system of identical fermions in thermodynamic equilibrium, the average number of fermions
## in a single-particle state `i` is given by the Fermi-Dirac distribution:

# Law: N_i = 1 / (exp((E_i - mu)/(k * T)) + 1)
## N_i - average number of fermions in single-particle state `i`,
##       also known as occupancy of energy state `i`
## E_i - energy of single-particle state `i`
## mu - total chemical potential
## k - Boltzmann constant
## T - absolute temperature

# Note
## If the energy states are degenerate, i.e. two or more particles are on the same energy level,
## the average number of fermions can be found by multiplying by the degeneracy `g_i` of the energy
## level.

occupancy_of_state = Symbol("occupancy_of_state", dimensionless)
energy_of_state = Symbol("energy_of_state", units.energy)
total_chemical_potential = Symbol("total_chemical_potential", units.energy)
temperature = symbols.thermodynamics.temperature

law = Eq(
    occupancy_of_state, 1 / (exp(
    (energy_of_state - total_chemical_potential) / (units.boltzmann_constant * temperature)) + 1))


def print_law() -> str:
    return print_expression(law)


@validate_input(
    energy_of_state_=energy_of_state,
    total_chemical_potential_=total_chemical_potential,
    temperature_=temperature,
)
@validate_output(occupancy_of_state)
def calculate_occupancy_of_state(
    energy_of_state_: Quantity,
    total_chemical_potential_: Quantity,
    temperature_: Quantity,
) -> float:
    result = law.rhs.subs({
        energy_of_state: energy_of_state_,
        total_chemical_potential: total_chemical_potential_,
        temperature: temperature_,
    })
    return convert_to_float(result)
