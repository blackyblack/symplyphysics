from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.optics import optical_strength_of_spherical_lens_from_refractive_indices_of_environment_and_lens_and_focal_distances as spherical_lens_law
from symplyphysics.laws.optics import lens_focus_from_object_and_image as focus_law

# Description
#

# Law: F = R / 2
# Where:
## F - focus distance
## R - curvature radius of mirror

focus_distance = Symbol("focus_distance", units.length)
curvature_radius = Symbol("curvature_radius", units.length)

law = Eq(focus_distance, curvature_radius / 2)

# refraction index of air equal 1. Mirror has refraction index equal -n, where n - refraction index of environment
spherical_lens_eq = spherical_lens_law.law.subs({
    spherical_lens_law.curvature_radius_lens: curvature_radius,
    spherical_lens_law.refraction_index_lens: -1,
    spherical_lens_law.refraction_index_environment: 1
})

# Paste distances from object and image in law of spherical lens
spherical_lens_equation = spherical_lens_eq.subs({
    spherical_lens_law.distance_to_object: focus_law.distance_to_object,
    spherical_lens_law.distance_to_image: focus_law.distance_to_image,
})

# Multiply both sides of the equation by (-1)
spherical_lens_equation = Eq(spherical_lens_equation.lhs * (-1), spherical_lens_equation.rhs * (-1))

focus_equation = focus_law.law.subs({
    focus_law.focus_distance: focus_distance
})

assert expr_equals(spherical_lens_equation.lhs, focus_equation.rhs)

focus_eq = Eq(focus_equation.lhs, spherical_lens_equation.rhs)
focus_value = solve(focus_eq, focus_distance, dict=True)[0][focus_distance]

assert expr_equals(focus_value, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(curvature_radius_=curvature_radius)
@validate_output(focus_distance)
def calculate_focus_distance(curvature_radius_: Quantity) -> Quantity:
    solved = solve(law, focus_distance, dict=True)[0][focus_distance]
    result_expr = solved.subs({
        curvature_radius: curvature_radius_
    })
    result = Quantity(result_expr)
    return result
