from sympy import (Eq, solve, sqrt)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
)

# Description
## Newton's rings are annular interference maxima and minima that appear around the point of contact between
## a convex lens and a plane—parallel plate when light passes through the lens and plate.
## The order of interference indicates the number of whole wavelengths. When the difference in the path of the
## beam reflected from the convex surface of the lens at the glass-air boundary and the beam reflected from the
## plate at the air-glass boundary is equal to an entire wave "n", then these two waves hit the observation
## point with the same faces and demonstrate interference.

## Law is: r = sqrt(k * R * L / n), where
## r - radius of Newton's ring,
## k - order of interference,
## R - radius of curvature of lens
## L - wavelength in medium surrounding lens and plate,
## n - refractive index of medium between lens and plate.

radius = Symbol("radius", units.length)

wavelength = Symbol("wavelength", units.length)
order_interference = Symbol("order_interference", dimensionless)
radius_curvature = Symbol("radius_curvature", units.length)
refractive_index = Symbol("refractive_index", dimensionless)

law = Eq(radius, sqrt(order_interference * radius_curvature * wavelength / refractive_index))


def print_law() -> str:
    return print_expression(law)


@validate_input(wavelength_=wavelength,
    order_interference_=order_interference,
    refractive_index_=refractive_index,
    radius_curvature_=radius_curvature)
@validate_output(radius)
def calculate_radius(wavelength_: Quantity, order_interference_: float, refractive_index_: float, radius_curvature_: Quantity) -> Quantity:
    check = order_interference_ * 2 % 2
    if not (isinstance(order_interference_, int) or check == 1 or check == 0):
        raise ValueError("order_interference_ must be an integer or an odd number divided by 2.")
    if order_interference_ <= 0:
        raise ValueError("order_interference_ must be greater than 0.") 

    result_expr = solve(law, radius, dict=True)[0][radius]
    result_expr = result_expr.subs({
        wavelength: wavelength_,
        order_interference: order_interference_,
        refractive_index: refractive_index_,
        radius_curvature: radius_curvature_,
    })
    return Quantity(result_expr)
