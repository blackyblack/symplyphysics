from sympy import Eq
from symplyphysics import (
    clone_symbol,
    symbols,
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

# Law: dG = - S * dT + V * dp + mu * dN
## G - Gibbs free energy
## S - entropy
## T - absolute temperature
## V - volume
## p - pressure
## mu - chemical potential
## N - number of particles in system
## Notation: d(x) - full differential of `x`

# Note
## - Temperature `T`, pressure `P`, and particle count `N` are so called natural variables of Gibbs energy
##   as a thermodynamic potential. Note that `N` is the only extensive natural variable, the other two are intensive.
## - For a system with more than one type of particles, the last term can be represented as a sum over all
##   types of particles, i.e. `Sum(mu_i * d(N_i), i)`.

# Conditions
## - The system is in thermal equilibrium with its surroundings
## - There is only one type of particles in the system

gibbs_energy_change = Symbol("gibbs_energy_change", units.energy)
entropy = Symbol("entropy", units.energy / units.temperature)
temperature_change = clone_symbol(symbols.thermodynamics.temperature, "temperature_change")
volume = Symbol("volume", units.volume)
pressure_change = Symbol("pressure_change", units.pressure)
chemical_potential = Symbol("chemical_potential", units.energy)
particle_count_change = Symbol("particle_count_change", dimensionless)

law = Eq(
    gibbs_energy_change,
    -1 * entropy * temperature_change + volume * pressure_change +
    chemical_potential * particle_count_change,
)


def print_law() -> str:
    return print_expression(law)


#pylint: disable=too-many-arguments
@validate_input(
    entropy_=entropy,
    temperature_change_=temperature_change,
    volume_=volume,
    pressure_change_=pressure_change,
    chemical_potential_=chemical_potential,
    particle_count_change_=particle_count_change,
)
@validate_output(gibbs_energy_change)
def calculate_gibbs_energy_change(
    entropy_: Quantity,
    temperature_change_: Quantity,
    volume_: Quantity,
    pressure_change_: Quantity,
    chemical_potential_: Quantity,
    particle_count_change_: int,
) -> Quantity:
    result = law.rhs.subs({
        entropy: entropy_,
        temperature_change: temperature_change_,
        volume: volume_,
        pressure_change: pressure_change_,
        chemical_potential: chemical_potential_,
        particle_count_change: particle_count_change_,
    })
    return Quantity(result)
