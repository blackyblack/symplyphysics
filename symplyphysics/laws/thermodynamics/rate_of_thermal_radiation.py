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
## One of the methods an object can exchange energy with its environment is via thermal radiation by
## emitting or absorbing energy in the form of electromagnetic waves.

# Law: P = sigma * epsilon * A * T**4
## P - rate of energy emission/absorption
## sigma - Stefan-Boltzmann constant
## epsilon - emissivity of object's surface (1 for idealized black body radiator)
## A - surface area of object
## T - in case of emission, temperature of object's surface;
##     in case of absorption, temperature of the environment

energy_radiation_rate = Symbol("energy_radiation_rate", units.power)
surface_emissivity = Symbol("surface_emissivity", dimensionless)
surface_area = Symbol("surface_area", units.area)
temperature_emission_absorption = clone_symbol(symbols.thermodynamics.temperature,
    "temperature_emission_absorption")

law = Eq(
    energy_radiation_rate,
    units.stefan_boltzmann_constant * surface_emissivity * surface_area *
    temperature_emission_absorption**4,
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    surface_emissivity_=surface_emissivity,
    surface_area_=surface_area,
    temperature_=temperature_emission_absorption,
)
@validate_output(energy_radiation_rate)
def calculate_energy_radiation_rate(
    surface_emissivity_: float,
    surface_area_: Quantity,
    temperature_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        surface_emissivity: surface_emissivity_,
        surface_area: surface_area_,
        temperature_emission_absorption: temperature_,
    })
    return Quantity(result)
