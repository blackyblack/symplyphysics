from sympy import Eq, solve
from symplyphysics import units, Quantity, Symbol, print_expression, validate_input, validate_output

# Description
# Law: F = ρ * g * V
# p - density of liquid
# g - acceleration due to gravity
# V - volume of the submerged part of the body

liquid_density = Symbol("liquid_density", units.mass / units.volume)
g = Symbol("g", units.acceleration)
submerged_volume = Symbol("submerged_volume", units.volume)
archimedes_force = Symbol("archimedes_force", units.force)

law = Eq(archimedes_force, liquid_density * g * submerged_volume)

def print_law() -> str:
    return print_expression(law)

@validate_input(density_=liquid_density, volume_=submerged_volume, g_=g)
@validate_output(archimedes_force)
def calculate_archimedes_force(density_: Quantity, volume_: Quantity, g_: Quantity) -> Quantity:
    result_force_expr = solve(law, archimedes_force, dict=True)[0][archimedes_force]
    result_expr = result_force_expr.subs({liquid_density: density_, submerged_volume: volume_, g: g_})
    return Quantity(result_expr)
