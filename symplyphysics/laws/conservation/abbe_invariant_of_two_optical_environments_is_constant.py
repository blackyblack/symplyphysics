"""
Abbe invariant of two optical environments is constant
======================================================

The point :math:`S` is located on the front of the optical axis,
i.e. on the part that is outside the spherical lens (outside).
The point :math:`S'` is located on a part of the optical axis inside the lens
The Abbe's invariant connects the front and back segments :math:`S` and :math:`S'`,
allowing one of them to be determined if the second one is known

**Conditions:**

#. Abbe's formula is valid only for paraxial rays;
#. Law is valid for one refractive surface
#. All rays emanating from point :math:`S` and forming different but necessarily small angles with the
   axis will pass through point :math:`S'` after refraction.

**Links:**

#. `OptoWiki <https://www.optowiki.info/glossary/abbe-invariant/>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

curvature_radius = symbols.radius_of_curvature
"""
:symbols:`radius_of_curvature`.
"""

medium_refraction_index = clone_as_symbol(symbols.relative_refractive_index, subscript="0")
"""
:symbols:`relative_refractive_index` of the medium.
"""

distance_from_object = clone_as_symbol(symbols.euclidean_distance, display_symbol="d_o", display_latex="d_\\text{o}")
"""
:symbols:`euclidean_distance` from lens to object.
"""

lens_refraction_index = symbols.relative_refractive_index
"""
:symbols:`relative_refractive_index` of the lens material.
"""

distance_from_image = clone_as_symbol(symbols.euclidean_distance, display_symbol="d_i", display_latex="d_\\text{i}")
"""
:symbols:`euclidean_distance` from lens to image.
"""

abbe_invariant_environment = medium_refraction_index * ((1 / distance_from_object) -
    (1 / curvature_radius))
abbe_invariant_lens = lens_refraction_index * ((1 / distance_from_image) - (1 / curvature_radius))

law = Eq(abbe_invariant_environment, abbe_invariant_lens)
"""
:laws:symbol::

:laws:latex::
"""

# NOTE:
## proofs: https://studme.org/341451/matematika_himiya_fizik/prelomlenie_otrazhenie_sveta_sfericheskoy_poverhnosti



@validate_input(
    curvature_radius_=curvature_radius,
    refraction_index_environment_=medium_refraction_index,
    distance_from_object_=distance_from_object,
    distance_from_image_=distance_from_image,
)
@validate_output(lens_refraction_index)
def calculate_refraction_index_lens(
    distance_from_object_: Quantity,
    distance_from_image_: Quantity,
    curvature_radius_: Quantity,
    refraction_index_environment_: float,
) -> float:
    solved = solve(law, lens_refraction_index, dict=True)[0][lens_refraction_index]
    result_expr = solved.subs({
        curvature_radius: curvature_radius_,
        medium_refraction_index: refraction_index_environment_,
        distance_from_object: distance_from_object_,
        distance_from_image: distance_from_image_,
    })
    return convert_to_float(result_expr)
