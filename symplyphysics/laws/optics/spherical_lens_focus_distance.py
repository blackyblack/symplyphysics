from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression,
    validate_input, validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.optics import lens_focus_from_object_and_image as focus_law

# Description
## According to the definition of the focus of an ellipse,
## the center of which is at the origin, the focus has coordinates F(c, 0),
## where c^2 = a^2 + b^2 ("a" and "b" are the main semi-axes of the ellipse).
## For the special case of an ellipse (circle): a = b = R, respectively,
## the focal length in projection on any of the axes of the coordinate system:

# Law: F = R/2
# Where:
## F - focus distance for spherical lens
## R - curvature radius of lens

focus_distance = Symbol("focus_distance", units.length)
radius_lens = Symbol("radius_lens", units.length)

law = Eq(focus_distance, radius_lens / 2)

# The focus of a spherical mirror is at a point from which the distance to the object
# is equal to the distance to the image - and is equal curvature radius
focus_eq = focus_law.law.subs({
    focus_law.distance_to_image: radius_lens,
    focus_law.distance_to_object: radius_lens,
    focus_law.focus_distance: focus_distance
})

focus_value = solve(focus_eq, focus_distance, dict=True)[0][focus_distance]
assert expr_equals(focus_value, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(radius_lens_=radius_lens)
@validate_output(focus_distance)
def calculate_focus_distance(radius_lens_: Quantity) -> Quantity:
    solved = solve(law, focus_distance, dict=True)[0][focus_distance]
    result_expr = solved.subs({
        radius_lens: radius_lens_
    })
    result = Quantity(result_expr)
    return result
