from sympy import (Eq, solve, sin)
from symplyphysics import (units, Quantity, Symbol, print_expression, dimensionless, validate_input,
    validate_output, angle_type)

# Description
## Diffraction from a three-dimensional periodic structure such as atoms in a crystal is called Bragg
## diffraction. This is similar to what happens when waves are scattered on a diffraction grating. Bragg
## diffraction is a consequence of interference between waves reflected from crystal planes.

## Law: d = n * L / (2 * sin(O)), where
## d - diffraction_order,
## n - refractive index of lens material,
## L - wavelength,
## O - sliding angle.

distance = Symbol("distance", units.length)

diffraction_order = Symbol("diffraction_order", dimensionless)
wavelength = Symbol("wavelength", units.length)
angle = Symbol("angle", angle_type)

law = Eq(distance, diffraction_order * wavelength / (2 * sin(angle)))


def print_law() -> str:
    return print_expression(law)


@validate_input(diffraction_order_=diffraction_order,
    wavelength_=wavelength,
    angle_=angle,)
@validate_output(distance)
def calculate_distance(diffraction_order_: float, wavelength_: Quantity, angle_: float | Quantity) -> Quantity:
    result_expr = solve(law, distance, dict=True)[0][distance]
    distance_applied = result_expr.subs({
        diffraction_order: diffraction_order_,
        wavelength: wavelength_,
        angle: angle_
    })
    return Quantity(distance_applied)
