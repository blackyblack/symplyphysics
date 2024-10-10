"""
Focal length of a concave spherical mirror
==========================================

The focal length of a spherical mirror is equal to half the radius of curvature.

**Notes:**

#. For a concave mirror, the focal length is positive.
#. For a convex mirror, the focal length is negative.
"""

from sympy import Eq, simplify, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.optics import (
    optical_strength_of_spherical_lens_from_refractive_indices_of_environment_and_lens_and_focal_distances as spherical_lens_law,
)
from symplyphysics.laws.optics import lens_focus_from_object_and_image as focus_law

focal_length = symbols.focal_length
"""
Mirror's :symbols:`focal_length`
"""

curvature_radius = symbols.radius_of_curvature
"""
Mirror's :symbols:`radius_of_curvature`
"""

law = Eq(focal_length, curvature_radius / 2)
"""
:laws:symbol::

:laws:latex::
"""

# Mirror has refraction index equal -n, where n - refraction index of environment
refraction_index = symbols.relative_refractive_index

spherical_lens_eq = spherical_lens_law.law.subs({
    spherical_lens_law.curvature_radius_lens: curvature_radius,
    spherical_lens_law.refraction_index_lens: -1 * refraction_index,
    spherical_lens_law.refraction_index_environment: refraction_index
})

# Paste distances from object and image in law of spherical lens
spherical_lens_equation = spherical_lens_eq.subs({
    spherical_lens_law.distance_to_object: focus_law.distance_to_object,
    spherical_lens_law.distance_to_image: focus_law.distance_to_image,
})

# Divide both sides of equation to refraction index
spherical_lens_equation = simplify(
    Eq(spherical_lens_equation.lhs / refraction_index,
    spherical_lens_equation.rhs / refraction_index))
spherical_lens_equation = Eq(curvature_radius / 2,
    solve(spherical_lens_equation, curvature_radius)[0] / 2)

focus_equation = focus_law.law.subs({focus_law.focus_distance: focal_length})
focus_equation = Eq(focal_length, solve(focus_equation, focal_length)[0])

focus_value = solve([focus_equation, spherical_lens_equation],
    (focal_length, focus_law.distance_to_image * focus_law.distance_to_object /
    (focus_law.distance_to_image + focus_law.distance_to_object)),
    dict=True)[0][focal_length]
assert expr_equals(focus_value, law.rhs)


@validate_input(curvature_radius_=curvature_radius)
@validate_output(focal_length)
def calculate_focus_distance(curvature_radius_: Quantity) -> Quantity:
    solved = solve(law, focal_length, dict=True)[0][focal_length]
    result_expr = solved.subs({curvature_radius: curvature_radius_})
    return Quantity(result_expr)
