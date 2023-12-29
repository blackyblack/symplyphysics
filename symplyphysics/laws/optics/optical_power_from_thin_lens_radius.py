from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, dimensionless, validate_input,
    validate_output)

# Description
## Optical power is the degree to which a lens, mirror, or other optical system converges or diverges light.
## For two or more thin lenses close together, the optical power of the combined lenses is approximately
## equal to the sum of the optical powers of each lens.

## Law: P = (n - n0)(1/R1 - 1/R2)
## Where:
## P - optical power of thin lense,
## n - refractive index of lens material,
## n0 - refractive index of medium,
## R1 - radius of the front surface of lens,
## R2 - radius of the back surface of lens.

# Conditions
## It is important to consider the sign of the radius of the surface!
## If the front surface of the lens is convex, then R1 > 0, and if concave, then R1 < 0.
## Conversely, for the back surface: if convex, then R2 < 0, if concave, then R2 > 0.
## Thus, for a biconvex lens, R1 > 0, R2 < 0, and for a biconcave lens, R1 < 0, R2 > 0.

optical_power = Symbol("optical_power", 1 / units.length)
lense_refractive_index = Symbol("lense_refractive_index", dimensionless)
medium_refractive_index = Symbol("medium_refractive_index", dimensionless)
front_radius = Symbol("front_radius", units.length)
back_radius = Symbol("back_radius", units.length)

law = Eq(optical_power,
    (lense_refractive_index - medium_refractive_index) * (1 / front_radius - 1 / back_radius))


def print_law() -> str:
    return print_expression(law)


@validate_input(object_lense_refractive_index_=lense_refractive_index,
    medium_refractive_index_=medium_refractive_index,
    front_radius_=front_radius,
    back_radius_=back_radius)
@validate_output(optical_power)
def calculate_optical_power(object_lense_refractive_index_: float, medium_refractive_index_: float,
    front_radius_: Quantity, back_radius_: Quantity) -> Quantity:
    result_expr = solve(law, optical_power, dict=True)[0][optical_power]
    optical_power_applied = result_expr.subs({
        lense_refractive_index: object_lense_refractive_index_,
        medium_refractive_index: medium_refractive_index_,
        front_radius: front_radius_,
        back_radius: back_radius_
    })
    return Quantity(optical_power_applied)