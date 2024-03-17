from sympy import (
    Eq,
    solve,
)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    angle_type,
)

# Description
## The resolution of the telescope is the minimum angular distance between point objects that can be
## distinguished separately in the telescope.

## Law is: r = 1.22 * (L / D), where
## r - resolution of telescope,
## L - wavelength,
## D - lens diameter.

resolution = Symbol("resolution", angle_type)

wavelength = Symbol("wavelength", units.length)
lens_diameter = Symbol("lens_diameter", units.length)

constant_rad = Quantity(1.22 * units.radian)

law = Eq(resolution, constant_rad * (wavelength / lens_diameter))


def print_law() -> str:
    return print_expression(law)


@validate_input(wavelength_=wavelength, lens_diameter_=lens_diameter)
@validate_output(resolution)
def calculate_resolution(wavelength_: Quantity, lens_diameter_: Quantity) -> Quantity:
    result_expr = solve(law, resolution, dict=True)[0][resolution]
    result_expr = result_expr.subs({
        wavelength: wavelength_,
        lens_diameter: lens_diameter_,
    })
    return Quantity(result_expr)
