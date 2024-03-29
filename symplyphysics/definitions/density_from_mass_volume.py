from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, symbols)

# Description
## The density (more precisely, the volumetric mass density), of a substance
## is its mass per unit volume.

# Definition: ρ = m / V
# Where:
## m is the mass
## V is volume
## ρ is the density

volume = Symbol("volume", units.volume)
density = Symbol("density", units.mass / units.volume)

definition = Eq(density, symbols.basic.mass / volume)

definition_units_SI = units.kilogram / units.meter**3


def print_law() -> str:
    return print_expression(definition)


@validate_input(mass_=symbols.basic.mass, volume_=volume)
@validate_output(density)
def calculate_density(mass_: Quantity, volume_: Quantity) -> Quantity:
    solved = solve(definition, density, dict=True)[0][density]
    result_expr = solved.subs({symbols.basic.mass: mass_, volume: volume_})
    return Quantity(result_expr)
