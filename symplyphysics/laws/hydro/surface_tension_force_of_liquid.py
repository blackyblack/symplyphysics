from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The surface tension force is directed tangentially to the surface of the liquid, perpendicular to
## the section of the contour on which it acts and is proportional to the length of this section.

## Law is: F = g * l, where
## F - surface tension force,
## g - coefficient of surface tension of the liquid,
## l - length of the liquid contour.

force = Symbol("force", units.force)

surface_coefficient = Symbol("surface_coefficient", units.force / units.length)
contour_length = Symbol("contour_length", units.length)

law = Eq(force, surface_coefficient * contour_length)


def print_law() -> str:
    return print_expression(law)


@validate_input(surface_coefficient_=surface_coefficient, contour_length_=contour_length)
@validate_output(force)
def calculate_force(surface_coefficient_: Quantity, contour_length_: Quantity) -> Quantity:
    result_expr = solve(law, force, dict=True)[0][force]
    result_expr = result_expr.subs({
        surface_coefficient: surface_coefficient_,
        contour_length: contour_length_,
    })
    return Quantity(result_expr)
