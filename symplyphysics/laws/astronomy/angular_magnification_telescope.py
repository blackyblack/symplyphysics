from sympy import (Eq, solve, S)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    convert_to,
    dimensionless,
)

# Description
## The angular magnification of the telescope shows how many times the angle at which an object is visible
## when viewed through a telescope is greater than when viewed with the eye.

## Law is: G = F / f, where
## G - angular magnification of telescope,
## F - focal length of the lens,
## f - focal length of the eyepiece.

angular_magnification = Symbol("angular_magnification", dimensionless)

focal_length_lens = Symbol("focal_length_lens", units.length)
focal_length_eyepiece = Symbol("focal_length_eyepiece", units.length)

law = Eq(angular_magnification, focal_length_lens / focal_length_eyepiece)


def print_law() -> str:
    return print_expression(law)


@validate_input(focal_length_lens_=focal_length_lens, focal_length_eyepiece_=focal_length_eyepiece)
@validate_output(angular_magnification)
def calculate_angular_magnification(focal_length_lens_: Quantity,
    focal_length_eyepiece_: Quantity) -> float:
    result_expr = solve(law, angular_magnification, dict=True)[0][angular_magnification]
    result_expr = result_expr.subs({
        focal_length_lens: focal_length_lens_,
        focal_length_eyepiece: focal_length_eyepiece_,
    })
    return float(convert_to(Quantity(result_expr), S.One).evalf())
