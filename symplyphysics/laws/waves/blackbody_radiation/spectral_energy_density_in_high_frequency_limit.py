from sympy import Eq, pi
from sympy.physics.units import speed_of_light, boltzmann_constant
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
## The Rayleigh-Jeans law is an approximation to the spectral radiance of electromagnetic radiation
## as a function of wave frequency from a blackbody at a given temperature through classical arguments.
## The Rayleigh-Jeans law agrees with experimental results at large wavelengths (i.e. at low frequencies)
## but strongly disagrees at short wavelengths (i.e. at high frequencies). This inconsistency is commonly
## known as ultraviolet catastrophe. See [Planck's law](./spectral_energy_density_at_all_frequencies.py) for
## the correct expression of radiation at all frequencies.

# Law: u_nu = 8 * pi * nu**2 * k * T / c**3
## u_nu - spectral energy density (energy per unit volume per unit frequency)
## nu - radiation frequency
## k - Boltzmann constant
## T - equilibrium temperature of black body
## c - speed of light

# Conditions
## - `h * nu << k * T`, i.e. the photon energy is much smaller than the thermal energy.

spectral_energy_density = Symbol("spectral_energy_density", units.energy / (units.volume * units.frequency))
radiation_frequency = Symbol("radiation_frequency", units.frequency)
equilibrium_temperature = clone_symbol(symbols.thermodynamics.temperature, "equilibrium_temperature")

law = Eq(
    spectral_energy_density,
    8 * pi * radiation_frequency**2 * boltzmann_constant * equilibrium_temperature / speed_of_light**3
)

# TODO: derive from Planck's law


def print_law() -> str:
    return print_expression(law)


@validate_input(
    radiation_frequency_=radiation_frequency,
    equilibrium_temperature_=equilibrium_temperature,
)
@validate_output(spectral_energy_density)
def calculate_spectral_energy_density(
    radiation_frequency_: Quantity,
    equilibrium_temperature_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        radiation_frequency: radiation_frequency_,
        equilibrium_temperature: equilibrium_temperature_,
    })
    return Quantity(result)
