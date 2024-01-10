from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
                           validate_output)

# Description
## The Archimedean force acting on a body immersed in a liquid (or gas) is equal to the weight of the liquid (or gas) displaced by the body.
## In this case, the weight of the body (i.e., the force with which the body acts on the support or suspension) immersed in the liquid decreases.

## Law: P_liquid = P * ( 1 - ⍴_liquid / ⍴)
## Where:
## P_liquid is body weight in liquid
## P is body weight in the air
## ⍴_liquid is the density of the liquid
## ⍴ is body density


weight_liquid = Symbol("weight_liquid", units.force)
weight_air = Symbol("weight_air", units.force)
liquid_density = Symbol("liquid_density", units.mass / units.volume)
body_density = Symbol("body_density", units.mass / units.volume)

law = Eq(weight_liquid, weight_air * (1 - (liquid_density / body_density)))


def print_law() -> str:
    return print_expression(law)


@validate_input(weight_air_=weight_air, liquid_density_=liquid_density, body_density_=body_density)
@validate_output(weight_liquid)
def calculate_weight(weight_air_: Quantity, liquid_density_, body_density_: Quantity) -> Quantity:
    result_expr = solve(law, weight_liquid, dict=True)[0][weight_liquid]
    result_weight = result_expr.subs({
        weight_air: weight_air_,
        liquid_density: liquid_density_,
        body_density: body_density_,
    })

    return Quantity(result_weight)
