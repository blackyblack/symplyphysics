from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, dimensionless, validate_input,
    validate_output)

# Description
## Optical power is the degree to which a lens, mirror, or other optical system converges or diverges light.
## For two or more thin lenses close together, the optical power of the combined lenses is approximately
## equal to the sum of the optical powers of each lens.

## Law: P = (n - n0)(1/R1 - 1/R2)
## Where:
## P - optical power of thin lens,
## n - refractive index of lens material,
## n0 - refractive index of medium,
## R1 - radius of the front surface of lens,
## R2 - radius of the back surface of lens.
## It is important to consider the sign of the radius of the surface!
## If the front surface of the lens is convex, then R1 > 0, and if concave, then R1 < 0.
## Conversely, for the back surface: if convex, then R2 < 0, if concave, then R2 > 0.
## Thus, for a biconvex lens, R1 > 0, R2 < 0, and for a biconcave lens, R1 < 0, R2 > 0.

# Conditions
## For front surface: R > 0 for convex, R < 0 for concave.
## For back surface: R < 0 for convex, R > 0 for concave.

# Links:
## Wikipedia, definition of optical power <https://en.wikipedia.org/wiki/Optical_power>
## Wikipedia, formula for effective focal length <https://en.wikipedia.org/wiki/Lens#Thin_lens_approximation>

optical_power = Symbol("optical_power", 1 / units.length)
lens_refractive_index = Symbol("lens_refractive_index", dimensionless)
medium_refractive_index = Symbol("medium_refractive_index", dimensionless)
front_radius = Symbol("front_radius", units.length)
back_radius = Symbol("back_radius", units.length)

law = Eq(optical_power,
    (lens_refractive_index - medium_refractive_index) * (1 / front_radius - 1 / back_radius))


@validate_input(lens_refractive_index_=lens_refractive_index,
    medium_refractive_index_=medium_refractive_index,
    front_radius_=front_radius,
    back_radius_=back_radius)
@validate_output(optical_power)
def calculate_optical_power(lens_refractive_index_: float, medium_refractive_index_: float,
    front_radius_: Quantity, back_radius_: Quantity) -> Quantity:
    result_expr = solve(law, optical_power, dict=True)[0][optical_power]
    optical_power_applied = result_expr.subs({
        lens_refractive_index: lens_refractive_index_,
        medium_refractive_index: medium_refractive_index_,
        front_radius: front_radius_,
        back_radius: back_radius_
    })
    return Quantity(optical_power_applied)
