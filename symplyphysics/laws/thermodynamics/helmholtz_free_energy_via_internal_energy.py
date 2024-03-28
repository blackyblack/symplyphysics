from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols,
)

# Description
## The Helmholtz free energy is a thermodynamic potential that measures the useful work obtainable 
## from a closed thermodynamic system at a constant temperature (isothermal).

# Law: F = U - T * S
## F - Helmholtz free energy
## U - internal energy of system
## T - absolute temperature of surroundings
## S - entropy of system

helmholtz_free_energy = Symbol("helmholtz_free_energy", units.energy)
internal_energy = Symbol("internal_energy", units.energy)
temperature = symbols.thermodynamics.temperature
entropy = Symbol("entropy", units.energy / units.temperature)

law = Eq(
    helmholtz_free_energy,
    internal_energy - temperature * entropy,
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    internal_energy_=internal_energy,
    temperature_=temperature,
    entropy_=entropy,
)
@validate_output(helmholtz_free_energy)
def calculate_helmholtz_free_energy(
    internal_energy_: Quantity,
    temperature_: Quantity,
    entropy_: Quantity,
) -> Quantity:
    # Note that internal energy and entropy are only known up to a constant

    result = law.rhs.subs({
        internal_energy: internal_energy_,
        temperature: temperature_,
        entropy: entropy_,
    })
    return Quantity(result)
