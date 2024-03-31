from sympy import (Eq, solve, S)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
    convert_to,
)

# Description
## The relative aperture of a telescope is the ratio of the diameter of the lens to its focal length. For visual observations,
## high-power telescopes give a larger exit pupil size, that is, the picture is bright and clear. A larger field of view allows
## you to observe extended objects, which include many galaxies and nebulae, that is, objects from Outer Space. In turn,
## non-high-power telescopes give a greater magnification, other things being equal, and are used in working with objects where
## details need to be considered, that is, with planets.

## Law is: A = D / F, where
## A - relative aperture of telescope,
## D - lens diameter,
## F - focal length of the lens.

relative_aperture = Symbol("relative_aperture", dimensionless)

lens_diameter = Symbol("lens_diameter", units.length)
focal_length_lens = Symbol("focal_length_lens", units.length)

law = Eq(relative_aperture, lens_diameter / focal_length_lens)


def print_law() -> str:
    return print_expression(law)


@validate_input(lens_diameter_=lens_diameter, focal_length_lens_=focal_length_lens)
@validate_output(relative_aperture)
def calculate_relative_aperture(lens_diameter_: Quantity, focal_length_lens_: Quantity) -> float:
    result_expr = solve(law, relative_aperture, dict=True)[0][relative_aperture]
    result_expr = result_expr.subs({
        lens_diameter: lens_diameter_,
        focal_length_lens: focal_length_lens_,
    })
    return float(convert_to(Quantity(result_expr), S.One).evalf())
