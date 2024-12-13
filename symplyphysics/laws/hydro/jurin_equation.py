from sympy import (Eq, solve, cos)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, angle_type)

# Description
## The Jurin equation determines the height of the liquid rise in the capillaries.
## In the absence of wetting  theta > 90, cos(theta) < 0, and the liquid level in the capillary drops by the value h. With full wetting  theta =0, cos (theta) = 1, and the radius of the meniscus is equal to the radius of the capillary.

## Conditions
## Surface of the meniscus is a sphere
## Height of the liquid is raised (lowered) h is much larger than the radius of the capillary r

## Law: h = 2 * a * cos(theta) / (⍴ * g * r)
## Where:
## h is the height of the liquid column
## a is the coefficient of surface tension of the liquid
## theta is the angle between the surface of a solid and the tangent to the surface of a liquid
## ⍴ is density of liquid
## g is acceleration of free fall
## r is capillary radius

# Links: Wikipedia <https://en.wikipedia.org/wiki/Jurin%27s_law>

height = Symbol("height", units.length)
surface_tension_coefficient = Symbol("surface_tension_coefficient", units.force / units.length)
angle = Symbol("angle", angle_type)
density_of_liquid = Symbol("density_of_liquid", units.mass / units.volume)
radius = Symbol("radius", units.length)

law = Eq(
    height, 2 * surface_tension_coefficient * cos(angle) /
    (density_of_liquid * radius * units.acceleration_due_to_gravity))


@validate_input(surface_tension_coefficient_=surface_tension_coefficient,
    angle_=angle,
    density_of_liquid_=density_of_liquid,
    radius_=radius)
@validate_output(height)
def calculate_height(surface_tension_coefficient_: Quantity, angle_: Quantity,
    density_of_liquid_: Quantity, radius_: Quantity) -> Quantity:
    result_expr = solve(law, height, dict=True)[0][height]
    result_height = result_expr.subs({
        surface_tension_coefficient: surface_tension_coefficient_,
        angle: angle_,
        density_of_liquid: density_of_liquid_,
        radius: radius_,
    })
    result_height_quantity = Quantity(result_height)
    if result_height_quantity.scale_factor < radius_.scale_factor:
        raise ValueError(
            f"The height must be greater than the radius. Currently {result_height_quantity.scale_factor} < {radius_.scale_factor}"
        )
    return result_height_quantity
