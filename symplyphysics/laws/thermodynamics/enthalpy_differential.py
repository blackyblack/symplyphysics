from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols,
)

# Description
## The fundamental thermodynamic relation are fundamental equations which demonstate how important
## thermodynamic quantities depend on variables that are measurable experimentally.

# Law: dH = T * dS + V * dp + mu * dN
## H - [enthalpy](./enthalpy_is_internal_energy_plus_pressure_energy.py)
## T - absolute temperature
## S - entropy
## p - pressure
## V - volume
## mu - chemical potential
## N - number of particles in system
## Notation: `d(x)` is an exact differential of `x`

# Note
## - Entropy `S`, pressure `p`, and particle count `N` are so called natural variables of enthalpy as a
##   thermodynamic potential
## - For a system with more than one type of particles, the last term can be represented as a sum over all
##   types of particles, i.e. `Sum(mu_i * d(N_i), i)`.

# Conditions
## - The system is in thermal equilibrium with the environment
## - Only reversible prossesses or pure heat transfer are considered
## - There is only one type of particles in the system.

enthalpy_change = Symbol("enthalpy_change", units.energy)
temperature = symbols.thermodynamics.temperature
entropy_change = Symbol("entropy_change", units.energy / units.temperature)
volume = Symbol("volume", units.volume)
pressure_change = Symbol("pressure_change", units.pressure)
chemical_potential = Symbol("chemical_potential", units.energy)
particle_count_change = Symbol("particle_count_change", dimensionless)

law = Eq(
    enthalpy_change,
    temperature * entropy_change + volume * pressure_change + chemical_potential * particle_count_change,
)

# TODO: Derive from definition of enthalpy and internal energy differential


def print_law() -> str:
    return print_expression(law)


@validate_input(
    temperature_=temperature,
    entropy_change_=entropy_change,
    volume_=volume,
    pressure_change_=pressure_change,
    chemical_potential_=chemical_potential,
    particle_count_change_=particle_count_change,
)
@validate_output(enthalpy_change)
# pylint: disable=too-many-arguments
def calculate_enthalpy_change(
    temperature_: Quantity,
    entropy_change_: Quantity,
    volume_: Quantity,
    pressure_change_: Quantity,
    chemical_potential_: Quantity,
    particle_count_change_: int,
) -> Quantity:
    result = law.rhs.subs({
        temperature: temperature_,
        entropy_change: entropy_change_,
        volume: volume_,
        pressure_change: pressure_change_,
        chemical_potential: chemical_potential_,
        particle_count_change: particle_count_change_,
    })
    return Quantity(result)
