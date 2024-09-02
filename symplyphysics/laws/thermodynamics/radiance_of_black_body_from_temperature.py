from sympy import (Eq, solve)
from symplyphysics import (symbols, units, Quantity, Symbol, print_expression, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import (
    radiant_exitance_is_radiant_flux_emitted_per_unit_area as exitance_def,)
from symplyphysics.laws.thermodynamics import radiation_power_via_temperature as radiation_law

# Description
## The Stefan–Boltzmann law, also known as Stefan's law, states that the total energy radiated per
## unit surface area of a black body in unit time (known variously as the black-body irradiance,
## energy flux density, radiant flux, or the emissive power), j*, is directly proportional to the fourth
## power of the black body's thermodynamic temperature T (also called absolute temperature).

# Law: j* = sigma*T^4, where
## j* is radiant heat energy (or radiant exitance, i.e. energy radiated per unit surface area per unit time),
## sigma is constant of proportionality, called the Stefan–Boltzmann constant,
## T is temperature of a completely black body

# Note
## j* = epsilon*sigma*T^4, where ε is the integral absorption capacity of the body. For a completely black body ε = 1.

radiance = Symbol("radiance", units.power / units.area)
temperature = symbols.thermodynamics.temperature

law = Eq(radiance, units.stefan_boltzmann_constant * temperature**4)

# Derive from law of thermal radiation power

_thermal_radiation_power = radiation_law.law.rhs.subs({
    radiation_law.emissivity: 1,  # see note, epsilon = 1 for idealized black body
    radiation_law.temperature: temperature,
    radiation_law.surface_area: exitance_def.area,
})

_radiant_exitance_derived = exitance_def.definition.rhs.subs(
    exitance_def.radiant_flux(exitance_def.area),
    _thermal_radiation_power,
).doit()

assert expr_equals(_radiant_exitance_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(temperature_=temperature)
@validate_output(radiance)
def calculate_radiance(temperature_: Quantity) -> Quantity:
    solved = solve(law, radiance, dict=True)[0][radiance]
    result_expr = solved.subs(temperature, temperature_)
    return Quantity(result_expr)
