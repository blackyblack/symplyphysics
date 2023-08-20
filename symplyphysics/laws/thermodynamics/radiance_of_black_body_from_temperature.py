from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## The Stefan–Boltzmann law, also known as Stefan's law, states that the total energy radiated per
## unit surface area of a black body in unit time (known variously as the black-body irradiance,
## energy flux density, radiant flux, or the emissive power), j*, is directly proportional to the fourth
## power of the black body's thermodynamic temperature T (also called absolute temperature).

# Law: j* = sigma*T^4, where
## j* is radiant heat energy,
## sigma is constant of proportionality, called the Stefan–Boltzmann constant,
## T is temperature of a completely black body

# Note
## j* = epsilon*sigma*T^4, where ε is the integral absorption capacity of the body. For a completely black body ε = 1.

radiance = Symbol("radiance", units.power / units.area)
temperature = Symbol("temperature", units.temperature)

law = Eq(radiance, units.stefan_boltzmann_constant * (temperature**4))


def print_law() -> str:
    return print_expression(law)


@validate_input(temperature_=temperature)
@validate_output(radiance)
def calculate_radiance(temperature_: Quantity) -> Quantity:
    solved = solve(law, radiance, dict=True)[0][radiance]
    result_expr = solved.subs(temperature, temperature_)
    return Quantity(result_expr)
