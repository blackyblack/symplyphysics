from sympy import Eq, exp, pi
from sympy.physics.units import planck, speed_of_light, boltzmann_constant
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
## Planck's radiation law describes the spectral density of electromagnetic radiation emitted
## by a black body in thermal equlibrium at a given temperature when there is no net flow of
## matter or energy between the body and its environment.

# Law: u_nu = (8 * pi * h * nu**3 / c**3) / (exp(h * nu / (k * T)) - 1))
## u_nu - spectral energy density (energy per unit volume per unit frequency)
## h - Planck constant
## nu - radiation frequency
## c - speed of light
## k - Boltzmann constant
## T - equilibrium temperature of black body

spectral_energy_density = Symbol("spectral_energy_density", units.energy / (units.volume * units.frequency))
radiation_frequency = Symbol("radiation_frequency", units.frequency)
equilibrium_temperature = clone_symbol(symbols.thermodynamics.temperature, "equilibrium_temperature")

law = Eq(
    spectral_energy_density,
    (8 * pi * planck * radiation_frequency**3 / speed_of_light**3)
    / (exp(planck * radiation_frequency / (boltzmann_constant * equilibrium_temperature)) - 1)
)


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
