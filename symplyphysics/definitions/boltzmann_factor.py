from sympy import Eq, exp, S
from symplyphysics import (
    convert_to,
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
## The Boltzmann factor is an exponential factor that appears in many formulas of statistical physics
## and thermodynamics, e.g. the canonical partition function of a classical discrete system.

# Definition: f = exp(-E_i / (k * T))
## f - Boltzmann factor
## E_i - energy of state i
## k - Boltzmann constant
## T - equilibrium temperature

boltzmann_factor = Symbol("boltzmann_factor", dimensionless)
energy_of_state = Symbol("energy_of_state", units.energy)
equilibrium_temperature = clone_symbol(symbols.thermodynamics.temperature,
    "equilibrium_temperature")

definition = Eq(boltzmann_factor,
    exp(-1 * energy_of_state / (units.boltzmann_constant * equilibrium_temperature)))


def print_law() -> str:
    return print_expression(definition)


@validate_input(
    energy_of_state_=energy_of_state,
    equilibrium_temperature_=equilibrium_temperature,
)
@validate_output(boltzmann_factor)
def calculate_boltzmann_factor(
    energy_of_state_: Quantity,
    equilibrium_temperature_: Quantity,
) -> float:
    result = definition.rhs.subs({
        energy_of_state: energy_of_state_,
        equilibrium_temperature: equilibrium_temperature_,
    })
    return float(convert_to(Quantity(result), S.One))
