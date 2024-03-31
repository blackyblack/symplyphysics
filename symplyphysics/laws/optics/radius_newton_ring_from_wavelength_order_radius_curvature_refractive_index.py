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
## r - radius of dark Newton's ring,
## k - order of ring,
## R - radius of curvature of lens
## L - wavelength in a vacuum,
## n - refractive index of medium between lens and plate.

radius = Symbol("radius", units.length)

wavelength = Symbol("wavelength", units.length)
order_of_ring = Symbol("order_of_ring", dimensionless)
radius_curvature = Symbol("radius_curvature", units.length)
refractive_index_between_lens_plate = Symbol("refractive_index_between_lens_plate", dimensionless)

law = Eq(radius, sqrt(order_of_ring * radius_curvature * wavelength / refractive_index_between_lens_plate))


def print_law() -> str:
    return print_expression(law)


@validate_input(wavelength_=wavelength,
    order_of_ring_=order_of_ring,
    refractive_index_between_lens_plate_=refractive_index_between_lens_plate,
    radius_curvature_=radius_curvature)
@validate_output(radius)
def calculate_radius(wavelength_: Quantity, order_of_ring_: int, refractive_index_between_lens_plate_: float, radius_curvature_: Quantity) -> Quantity:
    if order_of_ring_ <= 0:
        raise ValueError("Order of ring must be greater than 0.") 

    result_expr = solve(law, radius, dict=True)[0][radius]
    result_expr = result_expr.subs({
        wavelength: wavelength_,
        order_of_ring: order_of_ring_,
        refractive_index_between_lens_plate: refractive_index_between_lens_plate_,
        radius_curvature: radius_curvature_,
    })
    return Quantity(result_expr)
