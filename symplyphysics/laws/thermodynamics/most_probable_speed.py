from sympy import Eq, sqrt
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
## The most probable speed is the speed at which the Maxwell-Boltzmann speed distribution function
## is maximum.

# Law: v_prob = sqrt(2 * k * T / m)
## v_prob - most probable speed
## k - Boltzmann constant
## T - equilibrium temperature
## m - particle mass

most_probable_speed = Symbol("most_probable_speed", units.velocity, positive=True)
equilibrium_temperature = clone_symbol(symbols.thermodynamics.temperature, "equilibrium_temperature", positive=True)
particle_mass = clone_symbol(symbols.basic.mass, "particle_mass", positive=True)

law = Eq(
    most_probable_speed,
    sqrt(2 * units.boltzmann_constant * equilibrium_temperature / particle_mass)
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    equilibrium_temperature_=equilibrium_temperature,
    particle_mass_=particle_mass,
)
@validate_output(most_probable_speed)
def calculate_most_probable_speed(
    equilibrium_temperature_: Quantity,
    particle_mass_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        equilibrium_temperature: equilibrium_temperature_,
        particle_mass: particle_mass_,
    })
    return Quantity(result)
