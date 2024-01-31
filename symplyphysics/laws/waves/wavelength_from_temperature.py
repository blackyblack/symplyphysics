from sympy import Eq, solve
from sympy.physics.units import speed_of_light, planck, boltzmann
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input, validate_output)

# Description
## Any object with non-zero absolute temperature radiates energy. Absolutely black object is objects, wich doesn't reflect any radiation.
## Wavelength of most intesive radiation depends on object's temperature and is described by Win's law.

# Law: lambda = b/T, where
## lambda is wavelength,
## T is thermodinamical temperature of object,
## b is Win's constant, b = ch/ka, where c is speed of light, h is Planc's constant, k is Boltzman's constant and a = 4,965114.

intensive_wavelength = Symbol("intensive_wavelengh", units.length)
object_temperature = Symbol("object_temperature", units.temperature)
a = 4.965114
wins_constant = speed_of_light * planck / (boltzmann * a)

law = Eq(intensive_wavelength, wins_constant / object_temperature)


def print_law() -> str:
    return print_expression(law)

@validate_input(object_temperature_=object_temperature)
@validate_output(intensive_wavelength)
def calculate_intensive_wavelength(object_temperature_: Quantity) -> Quantity:
    result_wavelength_expr = solve(law, intensive_wavelength, dict=True)[0][intensive_wavelength]
    result_expr = result_wavelength_expr.subs({object_temperature: object_temperature_})
    return Quantity(result_expr)
