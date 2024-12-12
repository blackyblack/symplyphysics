from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
)

# Description
## The angular magnification of the telescope shows how many times the angle at which an object is visible
## when viewed through a telescope is greater than when viewed with the eye.

## Law is: G = F / f, where
## G - angular magnification of telescope,
## F - focal length of the lens,
## f - focal length of the eyepiece.

# Link: Physics LibreTexts, paragraph <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_III_-_Optics_and_Modern_Physics_(OpenStax)/02%3A_Geometric_Optics_and_Image_Formation/2.09%3A_Microscopes_and_Telescopes>

# Link: Physics LibreTexts, formula <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_III_-_Optics_and_Modern_Physics_(OpenStax)/02%3A_Geometric_Optics_and_Image_Formation/2.09%3A_Microscopes_and_Telescopes#mjx-eqn-eq2.36>

# TODO: move law to `optics`?

angular_magnification = Symbol("angular_magnification", dimensionless)

focal_length_lens = Symbol("focal_length_lens", units.length)
focal_length_eyepiece = Symbol("focal_length_eyepiece", units.length)

law = Eq(angular_magnification, focal_length_lens / focal_length_eyepiece)


@validate_input(focal_length_lens_=focal_length_lens, focal_length_eyepiece_=focal_length_eyepiece)
@validate_output(angular_magnification)
def calculate_angular_magnification(focal_length_lens_: Quantity,
    focal_length_eyepiece_: Quantity) -> float:
    result_expr = solve(law, angular_magnification, dict=True)[0][angular_magnification]
    result_expr = result_expr.subs({
        focal_length_lens: focal_length_lens_,
        focal_length_eyepiece: focal_length_eyepiece_,
    })
    return convert_to_float(result_expr)
