from sympy import (Eq, solve, S)
from symplyphysics import (units, Quantity, Symbol, print_expression,
    validate_input, validate_output, dimensionless, convert_to)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.laws.conservation import abbe_invariant_of_two_optical_environments_is_constant as abbe_conservation_law

# Description:
## The formula of the spherical refractive surface allows you to uniquely determine
## the position of the image if the position of the object is known and vice versa.

# Law: (n / (-lo)) + (n' / li) = (n' - n) / R
# Where:
## n - medium refraction index
## lo - distance from object to surface
## li - distance from image to surface
## R - surface curvature radius

# NOTES:
## - minus before lo if the origin of the coordinate system is placed at the intersection of the optical axis
##   and the surface of the lens with the radius lying on this optical axis;
## - proofs: https://studme.org/341451/matematika_himiya_fizik/prelomlenie_otrazhenie_sveta_sfericheskoy_poverhnosti.


distance_to_object = Symbol("distance_to_object", units.length)
distance_to_image = Symbol("distance_to_image", units.length)
curvature_radius_lens = Symbol("curvature_radius_lens", units.length)
refraction_index_environment = Symbol("refraction_index_environment", dimensionless)
refraction_index_lens = Symbol("refraction_index_lens", dimensionless)

law = Eq(
    (refraction_index_environment / (-distance_to_object)) + (refraction_index_lens / distance_to_image),
    (refraction_index_lens - refraction_index_environment) / curvature_radius_lens
)

# From Abbe's invariants:

abbe_invariant_environment = refraction_index_environment * ((1 / distance_to_object) - (1 / curvature_radius_lens))
abbe_invariant_lens = refraction_index_lens * ((1 / distance_to_image) - (1 / curvature_radius_lens))

invariant_conservation_eq = abbe_conservation_law.law.subs({
    abbe_conservation_law.abbe_invariant_environment: abbe_invariant_environment,
    abbe_conservation_law.abbe_invariant_lens: abbe_invariant_lens
})

radius_lens_from_law = solve(law, curvature_radius_lens, dict=True)[0][curvature_radius_lens]
radius_lens_from_invariants = solve(invariant_conservation_eq, curvature_radius_lens, dict=True)[0][curvature_radius_lens]
assert expr_equals(radius_lens_from_law, radius_lens_from_invariants)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    distance_to_object_=distance_to_object,
    distance_to_image_=distance_to_image,
    curvature_radius_lens_=curvature_radius_lens,
    refraction_index_environment_=refraction_index_environment
)
@validate_output(refraction_index_lens)
def calculate_refraction_index_lens(
        distance_to_object_: Quantity,
        distance_to_image_: Quantity,
        curvature_radius_lens_: Quantity,
        refraction_index_environment_: float
) -> float:
    solved = solve(law, refraction_index_lens, dict=True)[0][refraction_index_lens]
    result_expr = solved.subs({
        distance_to_object: distance_to_object_,
        distance_to_image: distance_to_image_,
        curvature_radius_lens: curvature_radius_lens_,
        refraction_index_environment: refraction_index_environment_
    })
    result = Quantity(result_expr)
    return float(convert_to(result, S.One).evalf())
