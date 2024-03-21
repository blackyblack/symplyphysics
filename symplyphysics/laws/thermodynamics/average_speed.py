from sympy import Eq, sqrt, pi
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols,
    clone_symbol,
)

# Description
## The average (mean) speed is the expected value of the speed distribution.

# Law: <v> = sqrt((8 / pi) * (k * T / m))
## <v> - average molecular speed
## k - Boltzmann constant
## T - equilibrium temperature
## m - molecular mass

# Conditions
## - Assuming the molecules are distributed according to Maxwell-Boltzmann statistics.

average_speed = Symbol("average_speed", units.velocity)
equilibrium_temperature = clone_symbol(symbols.thermodynamics.temperature, "equilibrium_temperature")
molecular_mass = clone_symbol(symbols.basic.mass, "molecular_mass")

law = Eq(
    average_speed, 
    sqrt(8 * units.boltzmann_constant * equilibrium_temperature / (pi * molecular_mass)),
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    equilibrium_temperature_=equilibrium_temperature,
    molecular_mass_=molecular_mass,
)
@validate_output(average_speed)
def calculate_average_speed(
    equilibrium_temperature_: Quantity,
    molecular_mass_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        equilibrium_temperature: equilibrium_temperature_,
        molecular_mass: molecular_mass_,
    })
    return Quantity(result)
