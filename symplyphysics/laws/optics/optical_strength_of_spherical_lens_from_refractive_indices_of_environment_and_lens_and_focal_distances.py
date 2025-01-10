"""
Optical power of spherical lens from refractive indeces and distances
=====================================================================

The formula of the spherical refractive surface allows you to uniquely determine
the position of the image if the position of the object is known and vice versa.

**Links:**

#. `Studme <https://studme.org/341451/matematika_himiya_fizik/prelomlenie_otrazhenie_sveta_sfericheskoy_poverhnosti>`__.

..
    TODO rename file
    TODO find English linj
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
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.conservation import abbe_invariant_of_two_optical_environments_is_constant as abbe_conservation_law

distance_to_object = clone_as_symbol(
    symbols.euclidean_distance,
    display_symbol="d_o",
    display_latex="d_\\text{o}",
)
"""
:symbols:`euclidean_distance` from lens to object.
"""

distance_to_image = clone_as_symbol(
    symbols.euclidean_distance,
    display_symbol="d_i",
    display_latex="d_\\text{i}",
)
"""
:symbols:`euclidean_distance` from lens to image.
"""

curvature_radius_lens = symbols.radius_of_curvature
"""
:symbols:`radius_of_curvature` of the lens surface.
"""

medium_refraction_index = clone_as_symbol(symbols.relative_refractive_index, subscript="0")
"""
:symbols:`relative_refractive_index` of the surrounding medium.
"""

lens_refraction_index = symbols.relative_refractive_index
"""
:symbols:`relative_refractive_index` of the lens material.
"""

law = Eq((medium_refraction_index /
    (-distance_to_object)) + (lens_refraction_index / distance_to_image),
    (lens_refraction_index - medium_refraction_index) / curvature_radius_lens)
"""
:laws:symbol::

:laws:latex::
"""

# TODO: prefix variables used in proof with underscore
# From Abbe's invariants:

_invariant_conservation_eq = abbe_conservation_law.law.subs({
    abbe_conservation_law.medium_refraction_index: medium_refraction_index,
    abbe_conservation_law.lens_refraction_index: lens_refraction_index,
    abbe_conservation_law.distance_from_image: distance_to_image,
    abbe_conservation_law.distance_from_object: distance_to_object,
    abbe_conservation_law.curvature_radius: curvature_radius_lens,
})

_radius_lens_from_law = solve(law, curvature_radius_lens, dict=True)[0][curvature_radius_lens]
_radius_lens_from_invariants = solve(_invariant_conservation_eq, curvature_radius_lens,
    dict=True)[0][curvature_radius_lens]
assert expr_equals(_radius_lens_from_law, _radius_lens_from_invariants)


@validate_input(distance_to_object_=distance_to_object,
    distance_to_image_=distance_to_image,
    curvature_radius_lens_=curvature_radius_lens,
    refraction_index_environment_=medium_refraction_index)
@validate_output(lens_refraction_index)
def calculate_refraction_index_lens(distance_to_object_: Quantity, distance_to_image_: Quantity,
    curvature_radius_lens_: Quantity, refraction_index_environment_: float) -> float:
    solved = solve(law, lens_refraction_index, dict=True)[0][lens_refraction_index]
    result_expr = solved.subs({
        distance_to_object: distance_to_object_,
        distance_to_image: distance_to_image_,
        curvature_radius_lens: curvature_radius_lens_,
        medium_refraction_index: refraction_index_environment_
    })
    return convert_to_float(result_expr)
