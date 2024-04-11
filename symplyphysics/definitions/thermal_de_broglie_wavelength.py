from sympy import Eq, sqrt, pi
from symplyphysics import (
    dimensionless,
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
## The thermal de Broglie wavelength can be roughly described as the average de Broglie wavelength
## of particles in an ideal gas at a specified temperature. When compared to average interparticle
## spacing in the gas, it can be used to tell if the gas can be considered to be a classical or
## Maxwell-Boltzmann gas, in which case the thermal wavelength must be much smaller than the average
## interparticle spacing. Otherwise, quantum effects must be taken into account.

# Definition: lambda = hbar * sqrt(2 * pi / (m * k * T))
## lambda - thermal de Broglie wavelength
## hbar - reduced Planck constant
## m - mass of gas particle
## k - Boltzmann constant
## T - gas temperature

thermal_wavelength = Symbol("thermal_wavelength", units.length)
particle_mass = clone_symbol(symbols.basic.mass, "particle_mass")
temperature = symbols.thermodynamics.temperature

law = Eq(
    thermal_wavelength,
    units.hbar * sqrt(2 * pi / (particle_mass * units.boltzmann_constant * temperature))
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    particle_mass_=particle_mass,
    temperature_=temperature,
)
@validate_output(thermal_wavelength)
def calculate_thermal_wavelength(
    particle_mass_: Quantity,
    temperature_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        particle_mass: particle_mass_,
        temperature: temperature_,
    })
    return Quantity(result)
