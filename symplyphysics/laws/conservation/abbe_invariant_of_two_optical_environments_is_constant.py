from sympy import Eq, solve
from symplyphysics import (
    Symbol,
    units,
    print_expression,
    Quantity,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
)

# Description
## The point S is located on the front of the optical axis,
## i.e. on the part that is outside the spherical lens (outside).
## The point S' is located on a part of the optical axis inside the lens
## The Abbe's invariant connects the front and back segments S and S',
## allowing one of them to be determined if the second one is known

# Law: n * ((1 / a) - (1 / R)) = n' * ((1 / b) - (1 / R))
# Where:
## n - refraction index in 1st optical environment
## n' - refraction index in 2nd optical environment
## a - distance from object to surface
## b - distance from image to surface
## R - surface curvature radius (the surface is convex)

# Conditions:
## - Abbe's formula is valid only for paraxial rays;
#№ - Law is valid for one refractive surface
## - All rays emanating from point S and forming different but necessarily small angles with the axis will pass through point S' after refraction.

# NOTE:
## proofs: https://studme.org/341451/matematika_himiya_fizik/prelomlenie_otrazhenie_sveta_sfericheskoy_poverhnosti

curvature_radius = Symbol("curvature_radius", units.length)

refraction_index_environment = Symbol("refraction_index_environment", dimensionless)
distance_from_object = Symbol("distance_from_object", units.length)

refraction_index_lens = Symbol("refraction_index_lens", dimensionless)
distance_from_image = Symbol("distance_from_image", units.length)

abbe_invariant_environment = refraction_index_environment * ((1 / distance_from_object) -
    (1 / curvature_radius))
abbe_invariant_lens = refraction_index_lens * ((1 / distance_from_image) - (1 / curvature_radius))

law = Eq(abbe_invariant_environment, abbe_invariant_lens)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    curvature_radius_=curvature_radius,
    refraction_index_environment_=refraction_index_environment,
    distance_from_object_=distance_from_object,
    distance_from_image_=distance_from_image,
)
@validate_output(refraction_index_lens)
def calculate_refraction_index_lens(
    distance_from_object_: Quantity,
    distance_from_image_: Quantity,
    curvature_radius_: Quantity,
    refraction_index_environment_: float,
) -> float:
    solved = solve(law, refraction_index_lens, dict=True)[0][refraction_index_lens]
    result_expr = solved.subs({
        curvature_radius: curvature_radius_,
        refraction_index_environment: refraction_index_environment_,
        distance_from_object: distance_from_object_,
        distance_from_image: distance_from_image_,
    })
    return convert_to_float(result_expr)
