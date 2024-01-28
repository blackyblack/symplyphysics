from sympy import (Eq, solve, S)
from symplyphysics import (units, Quantity, Symbol, print_expression,
    validate_input, validate_output, dimensionless, convert_to)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.laws.conservation import abbe_invariant_of_two_optical_environments_is_constant as abbe_conservation_law
from symplyphysics.laws.optics import abbe_invariant as abbe_invariant_law

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
lens_radius = Symbol("lens_radius", units.length)
refraction_index_medium = Symbol("refraction_index_medium", dimensionless)
refraction_index_lens = Symbol("refraction_index_lens", dimensionless)

law = Eq(
    (refraction_index_medium / (-distance_to_object)) + (refraction_index_lens / distance_to_image),
    (refraction_index_lens - refraction_index_medium) / lens_radius
)

# From Abbe's invariants:

medium_invariant_value = abbe_invariant_law.law.subs({
    abbe_invariant_law.curvature_radius: lens_radius,
    abbe_invariant_law.refraction_index: refraction_index_medium,
    abbe_invariant_law.distance_to_surface: distance_to_object
}).rhs
lens_invariant_value = abbe_invariant_law.law.subs({
    abbe_invariant_law.curvature_radius: lens_radius,
    abbe_invariant_law.refraction_index: refraction_index_lens,
    abbe_invariant_law.distance_to_surface: distance_to_image
}).rhs

invariant_conservation_eq = abbe_conservation_law.law.subs({
    abbe_conservation_law.invariant_abbe_before: medium_invariant_value,
    abbe_conservation_law.invariant_abbe_after: lens_invariant_value
})

radius_lens_from_law = solve(law, lens_radius, dict=True)[0][lens_radius]
radius_lens_from_invariants = solve(invariant_conservation_eq, lens_radius, dict=True)[0][lens_radius]
assert expr_equals(radius_lens_from_law, radius_lens_from_invariants)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    distance_to_object_=distance_to_object, distance_to_image_=distance_to_image,
    lens_radius_=lens_radius, refraction_index_medium_=refraction_index_medium
)
@validate_output(refraction_index_lens)
def calculate_refraction_index_lens(
        distance_to_object_: Quantity, distance_to_image_: Quantity,
        lens_radius_: Quantity, refraction_index_medium_: Quantity
) -> float:
    solved = solve(law, refraction_index_lens, dict=True)[0][refraction_index_lens]
    result_expr = solved.subs({
        distance_to_object: distance_to_object_,
        distance_to_image: distance_to_image_,
        lens_radius: lens_radius_,
        refraction_index_medium: refraction_index_medium_
    })
    result = Quantity(result_expr)
    return float(convert_to(result, S.One).evalf())
