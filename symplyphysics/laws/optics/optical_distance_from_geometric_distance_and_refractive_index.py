from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, dimensionless, validate_input,
    validate_output)

# Description
## Optical distance L is a geometric distance l in optical environment with refractive index n

# Law: L = n * l
# Where:
## L - optical distance
## n - refractive index of environment,
## l - geometric distance

optical_distance = Symbol("optical_power", units.length)
distance = Symbol("distance", units.length)
refractive_index = Symbol("refractive_index", dimensionless)

law = Eq(optical_distance, refractive_index * distance)


def print_law() -> str:
    return print_expression(law)


@validate_input(distance_=distance, refractive_index_=refractive_index)
@validate_output(optical_distance)
def calculate_optical_distance(distance_: float, refractive_index_: float) -> Quantity:
    result_expr = solve(law, optical_distance, dict=True)[0][optical_distance]
    optical_power_applied = result_expr.subs({
        distance: distance_,
        refractive_index: refractive_index_
    })
    return Quantity(optical_power_applied)
