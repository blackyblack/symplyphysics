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
## The formula of the internal energy differential can be integrated using the Euler's theorem on
## homogeneous functions to get the following expression:

# Law: U = T * S - p * V + mu * N
## U - internal energy
## T - temperature
## S - entropy
## p - pressure
## V - volume
## mu - chemical potential
## N - particle count

internal_energy = Symbol("internal_energy", units.energy)
temperature = symbols.thermodynamics.temperature
entropy = Symbol("entropy", units.energy / units.temperature)
pressure = Symbol("pressure", units.pressure)
volume = Symbol("volume", units.volume)
chemical_potential = Symbol("chemical_potential", units.energy)
particle_count = Symbol("particle_count", dimensionless)

law = Eq(
    internal_energy,
    temperature * entropy - pressure * volume + chemical_potential * particle_count
)

# TODO: derive from Euler's theorem and internal energy differential


def print_law() -> str:
    return print_expression(law)


@validate_input(
    temperature_=temperature,
    entropy_=entropy,
    pressure_=pressure,
    volume_=volume,
    chemical_potential_=chemical_potential,
    particle_count_=particle_count,
)
@validate_output(internal_energy)
# pylint: disable=too-many-arguments
def calculate_internal_energy(
    temperature_: Quantity,
    entropy_: Quantity,
    pressure_: Quantity,
    volume_: Quantity,
    chemical_potential_: Quantity,
    particle_count_: int,
) -> Quantity:
    result = law.rhs.subs({
        temperature: temperature_,
        entropy: entropy_,
        pressure: pressure_,
        volume: volume_,
        chemical_potential: chemical_potential_,
        particle_count: particle_count_,
    })
    return Quantity(result)
