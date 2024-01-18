from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
                           validate_output)

# Description
## If the body is on the surface of a liquid (floating), then only two forces act on it (Archimedes up and gravity down), which balance each other.

## Law: V_sub / V = ⍴ / ⍴_liquid
## Where:
## V_sub is the volume of the submerged body part
## V is full body volume
## ⍴ is body density
## ⍴_liquid is the density of the liquid

## Conditions
## The body is immersed in only one liquid.
## Gravity should be equal to the force of Archimedes.


submerged_volume = Symbol("submerged_volume", units.volume)
body_volume = Symbol("body_volume", units.volume)
body_density = Symbol("body_density", units.mass / units.volume)
liquid_density = Symbol("liquid_density", units.mass / units.volume)


law = Eq(submerged_volume / body_volume, body_density / liquid_density)


def print_law() -> str:
    return print_expression(law)


@validate_input(body_volume_=body_volume, body_density_=body_density, liquid_density_=liquid_density)
@validate_output(submerged_volume)
def calculate_submerged_volume(body_volume_: Quantity, body_density_, liquid_density_: Quantity) -> Quantity:
    result_expr = solve(law, submerged_volume, dict=True)[0][submerged_volume]
    result_volume = result_expr.subs({
        body_volume: body_volume_,
        body_density: body_density_,
        liquid_density: liquid_density_,
    })

    return Quantity(result_volume)
