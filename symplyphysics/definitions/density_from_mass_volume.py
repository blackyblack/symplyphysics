from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## The density (more precisely, the volumetric mass density), of a substance
## is its mass per unit volume.

## Definition: ρ = m / V
## Where:
## m is the mass
## V is volume
## ρ is the density

mass, volume, density = symbols('mass volume density')
definition = Eq(density, mass / volume)

definition_dimension_SI = units.kilogram / units.meter**3

def print():
    return pretty(definition, use_unicode=False)

def print_dimension():
    return pretty(definition_dimension_SI, use_unicode=False)

@validate_input(mass_=units.mass, volume_=units.volume)
@validate_output(units.mass / units.volume)
def calculate_density(mass_: Quantity, volume_: Quantity) -> Quantity:
    solved = solve(definition, density, dict=True)[0][density]
    result_expr = solved.subs({
        mass: mass_,
        volume: volume_})
    return expr_to_quantity(result_expr, 'volumetric_density')