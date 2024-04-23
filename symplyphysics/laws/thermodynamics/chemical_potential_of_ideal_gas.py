from sympy import Eq, log
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
## The chemical potential of an ideal gas can be calculated from its temperature, concentration,
## and thermal wavelength.

# Law: mu = k * T * log(n * lambda**3)
## mu - chemical potential of ideal gas
## k - Boltzmann constant
## T - gas temperature
## n - concentration of gas particles (particle count per unit volume)
## lambda - [thermal de Broglie wavelength](../../definitions/thermal_de_broglie_wavelength.py) of gas

chemical_potential = Symbol("chemical_potential", units.energy)
temperature = symbols.thermodynamics.temperature
concentration = Symbol("concentration", 1 / units.volume)
thermal_wavelength = Symbol("thermal_wavelength", units.length)

law = Eq(chemical_potential,
    units.boltzmann_constant * temperature * log(concentration * thermal_wavelength**3))


def print_law() -> str:
    return print_expression(law)


@validate_input(
    temperature_=temperature,
    concentration_=concentration,
    thermal_wavelength_=thermal_wavelength,
)
@validate_output(chemical_potential)
def calculate_chemical_potential(
    temperature_: Quantity,
    concentration_: Quantity,
    thermal_wavelength_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        temperature: temperature_,
        concentration: concentration_,
        thermal_wavelength: thermal_wavelength_,
    })
    return Quantity(result)
