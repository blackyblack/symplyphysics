from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, dimensionless, print_expression, validate_input,
    validate_output)

# Description
## The lens maker’s equation is a formula that gives us a relationship between
## the focal length, refractive index, and radii of curvature of the two spheres used in lenses.

## Law: 1 / F = (n - 1)(1/R1 - 1/R2)
## Where:
## F is the focal length (half the radius of curvature)
## n is the refractive index of the material used
## R1 is the radius of curvature of sphere 1
## R2 is the radius of curvature of sphere 2

# Conditions
## Lenses where the thickness is lesser, such that they are considered negligible
## in comparison to the radius of curvature, are referred to as thin lenses.=


focal_length = Symbol("focal_length", units.length)
refractive_index = Symbol("refractive_index", dimensionless)
radius_curvature_1 = Symbol("radius_curvature_1", units.length)
radius_curvature_2 = Symbol("radius_curvature_2", units.length)

law = Eq((1 / focal_length), (refractive_index - 1) * (1 / radius_curvature_1) - (1 / radius_curvature_2))


def print_law() -> str:
    return print_expression(law)


@validate_input(refractive_index_=refractive_index, radius_curvature_1_=radius_curvature_1,
                radius_curvature_2_=radius_curvature_2)
@validate_output(focal_length)
def calculate_length(refractive_index_: Quantity, radius_curvature_1_: Quantity, radius_curvature_2_: Quantity) -> Quantity:
    result_expr = solve(law, focal_length, dict=True)[0][focal_length]
    result_length = result_expr.subs({
        refractive_index: refractive_index_,
        radius_curvature_1: radius_curvature_1_,
        radius_curvature_2: radius_curvature_2_
    })

    return Quantity(result_length)