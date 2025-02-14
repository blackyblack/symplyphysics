from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
    symbols,
)

# Description
## The relative aperture of a telescope is the ratio of the diameter of the lens to its focal length. For visual observations,
## high-power telescopes give a larger exit pupil size, that is, the picture is bright and clear. A larger field of view allows
## you to observe extended objects, which include many galaxies and nebulae, that is, objects from Outer Space. In turn,
## non-high-power telescopes give a greater magnification, other things being equal, and are used in working with objects where
## details need to be considered, that is, with planets.

# Link: Wikipedia <https://en.wikipedia.org/wiki/F-number#Notation>

# TODO: move law to `optics`?
# TODO: update documentation

## Law is: A = D / F, where
## A - relative aperture of telescope,
## D - lens diameter,
## F - focal length of the lens.

relative_aperture = SymbolNew("A", dimensionless)

lens_diameter = symbols.diameter
focal_length_lens = symbols.focal_length

law = Eq(relative_aperture, lens_diameter / focal_length_lens)


@validate_input(lens_diameter_=lens_diameter, focal_length_lens_=focal_length_lens)
@validate_output(relative_aperture)
def calculate_relative_aperture(lens_diameter_: Quantity, focal_length_lens_: Quantity) -> float:
    result_expr = solve(law, relative_aperture, dict=True)[0][relative_aperture]
    result_expr = result_expr.subs({
        lens_diameter: lens_diameter_,
        focal_length_lens: focal_length_lens_,
    })
    return convert_to_float(result_expr)
