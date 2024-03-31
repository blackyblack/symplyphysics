from sympy import (Eq, solve, cos)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
    angle_type
)

# Description
## Interference in thin films is a phenomenon that occurs as a result of the separation of a light beam when
## reflected from the upper and lower boundaries of a thin film. As a result, there are two light waves that
## can interfere.
## The order of interference indicates the number of whole wavelengths. When the difference in the path of the
## beam reflected from the upper boundary of the film and the beam reflected from the lower boundary of the
## film is equal to an integer wave "n", then these two waves enter the observation point with the same phases
## and demonstrate interference.
## Order of interference can be chosen arbitrarily, to achieve the desired thin film thickness.
## Angle of refraction is angle between refracted ray and the normal.

## Law is: h = k * L / (2 * n * cos(O)), where
## h - film thickness,
## k - order of interference,
## L - wavelength in vacuum,
## n - refractive index of film
## O - angle of refraction.

film_thickness = Symbol("film_thickness", units.length)

wavelength = Symbol("wavelength", units.length)
order_interference = Symbol("order_interference", dimensionless)
angle_refraction = Symbol("angle_refraction", angle_type)
refractive_index = Symbol("refractive_index", dimensionless)

law = Eq(film_thickness, order_interference * wavelength / (2 * refractive_index * cos(angle_refraction)))


def print_law() -> str:
    return print_expression(law)


@validate_input(wavelength_=wavelength,
    order_interference_=order_interference,
    refractive_index_=refractive_index,
    angle_refraction_=angle_refraction)
@validate_output(film_thickness)
def calculate_film_thickness(wavelength_: Quantity, order_interference_: int, refractive_index_: float, angle_refraction_: float | Quantity) -> Quantity:
    if order_interference_ <= 0:
        raise ValueError("Order interference must be greater than 0.") 

    result_expr = solve(law, film_thickness, dict=True)[0][film_thickness]
    result_expr = result_expr.subs({
        wavelength: wavelength_,
        order_interference: order_interference_,
        refractive_index: refractive_index_,
        angle_refraction: angle_refraction_,
    })
    return Quantity(result_expr)
