from sympy import Eq, exp, S
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols,
    convert_to,
)

# Description
## Bose-Einstein statistics is a type of quantum statistics that applies to the physics of a
## system consisting of many non-interacting, identical particles that strictly do not obey
## the Pauli exclusion principle. Particles following Bose-Einstein statistics are called bosons
## and have integer values of spin.

# Law: n_i = 1 / (exp((E_i - mu)/(k * T)) - 1)
## n_i - occupancy of single-particle state `i`
## E_i - energy of single-particle state `i`
## mu - total chemical potential of system
## k - Boltzmann constant
## T - absolute temperature of the system

occupancy_of_state = Symbol("occupancy_of_state", dimensionless)
energy_of_state = Symbol("energy_of_state", units.energy)
total_chemical_potential = Symbol("total_chemical_potential", units.energy)
temperature = symbols.thermodynamics.temperature

law = Eq(
    occupancy_of_state, 1 / (exp(
    (energy_of_state - total_chemical_potential) / (units.boltzmann * temperature)) - 1))


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
    return float(convert_to(Quantity(result), S.One))
