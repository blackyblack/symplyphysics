from sympy import (Eq, solve, pi)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## Luminosity is an absolute measure of radiated electromagnetic energy (light) per unit time, and is synonymous with the radiant power emitted by a light-emitting object.
## In astronomy, luminosity is the total amount of electromagnetic energy emitted per unit of time by a star, galaxy, or other astronomical objects.

# Law: L = 4 * pi * R^2 * sigma * T^4

## L is luminosity
## sigma is constant of proportionality, called the Stefan–Boltzmann constant,
## R is the radius of the star
## T is the temperature of the star's photosphere

luminosity = Symbol("luminosity", units.power)
radius = Symbol("radius", units.length)
temperature = Symbol("temperature", units.temperature)

law = Eq(luminosity, 4 * pi * (radius**2) * units.stefan_boltzmann_constant * (temperature**4))


def print_law() -> str:
    return print_expression(law)


@validate_input(radius_=radius, temperature_=temperature)
@validate_output(luminosity)
def calculate_luminosity(radius_: Quantity, temperature_: Quantity) -> Quantity:
    solved = solve(law, luminosity, dict=True)[0][luminosity]
    result_expr = solved.subs({
        radius: radius_,
        temperature: temperature_
    })
    return Quantity(result_expr)
