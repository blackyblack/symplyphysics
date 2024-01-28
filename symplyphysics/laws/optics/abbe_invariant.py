from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression,
    validate_input, validate_output, dimensionless)

# Description
## The Abbe invariant is an optical characteristic of the refractive curve of the surface

# Law: Q = n * ((1 / d) - (1 / R))
# Where:
## n - surface refraction index
## d - distance to surface
## R - surface curvature radius
## Q - Abbe's invariant

# NOTE:
## proofs:
## https://studopedia.ru/10_198746_otritsatelnie-tendentsii-v-rabote-zao-tehnolog.html

refraction_index = Symbol("refraction_index", dimensionless)
curvature_radius = Symbol("curvature_radius", units.length)
distance_to_surface = Symbol("distance_to_surface", units.length)
abbe_invariant = Symbol("abbe_invariant", 1 / units.length)

law = Eq(abbe_invariant, refraction_index * ((1 / distance_to_surface) - (1 / curvature_radius)))


def print_law() -> str:
    return print_expression(law)


@validate_input(distance_to_surface_=distance_to_surface, curvature_radius_=curvature_radius, refraction_index_=refraction_index)
@validate_output(abbe_invariant)
def calculate_abbe_invariant(distance_to_surface_: Quantity, curvature_radius_: Quantity, refraction_index_: Quantity) -> Quantity:
    solved = solve(law, abbe_invariant, dict=True)[0][abbe_invariant]
    result_expr = solved.subs({
        refraction_index: refraction_index_,
        curvature_radius: curvature_radius_,
        distance_to_surface: distance_to_surface_
    })
    result = Quantity(result_expr)
    return result
