from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The fundamental thermodynamic relation are fundamental equations which demonstate how important
## thermodynamic quantities depend on variables that are measurable experimentally.

# Law: dF = -S * dT - p * dV + mu * dN
## F - Helmholtz free energy
## S - entropy
## T - absolute temperature
## p - pressure
## V - volume
## mu - chemical potential
## N - number of particles in the system
## Notation: d(x) - full differential of `x`

# Note
## - Temperature `T`, volume `V`, and particle count N are so called natural variables of Helmholtz free energy
##   as a thermodynamic potential.
## - For a system with more than one type of particles, the last term can be represented as a sum over all
##   types of particles, i.e. `Sum(mu_i * d(N_i), i)`.

# Conditions
## - The system is in thermal equilibrium with its surroundings
## - There is only one type of particles in the system

free_energy_change = Symbol("free_energy_change", units.energy)
entropy = Symbol("entropy", units.energy / units.temperature)
temperature_change = Symbol("temperature_change", units.temperature)
pressure = Symbol("pressure", units.pressure)
volume_change = Symbol("volume_change", units.volume)
chemical_potential = Symbol("chemical_potential", units.energy)
particle_count_change = Symbol("particle_count_change", dimensionless)

law = Eq(
    free_energy_change,
    -1 * entropy * temperature_change - pressure * volume_change + chemical_potential * particle_count_change,
)

# TODO: derive from the definition of Helmholtz free energy and internal energy differential


def print_law() -> str:
    return print_expression(law)


@validate_input(
    entropy_=entropy,
    temperature_change_=temperature_change,
    pressure_=pressure,
    volume_change_=volume_change,
    chemical_potential_=chemical_potential,
    particle_count_change_=particle_count_change,
)
@validate_output(free_energy_change)
def calculate_free_energy_change(
    entropy_: Quantity,
    temperature_change_: Quantity,
    pressure_: Quantity,
    volume_change_: Quantity,
    chemical_potential_: Quantity,
    particle_count_change_: int,
) -> Quantity:
    result = law.rhs.subs({
        entropy: entropy_,
        temperature_change: temperature_change_,
        pressure: pressure_,
        volume_change: volume_change_,
        chemical_potential: chemical_potential_,
        particle_count_change: particle_count_change_,
    })
    return Quantity(result)
