from sympy import pi
from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)
from symplyphysics.laws.nuclear.buckling import neutron_flux_for_uniform_parallelepiped as parallelepiped_flux

# Description
## Geometric buckling for the uniform parallelepiped reactor of physical dimensions a * b * c.
## The parallelepiped reactor is situated in cartesian geometry at the origin of coordinates.

## Law: Bg^2parallelepiped = (Pi / a)^2 + (Pi / b)^2 + (Pi / c)^2
## Where:
## Pi - Pi constant.
## a - width of the parallelepiped.
## b - length of the parallelepiped.
## c - height of the parallelepiped.
## Bg^2parallelepiped - squared geometric buckling for parallelepiped.
##   See [geometric buckling](./geometric_buckling_from_neutron_flux.py) implementation.

parallelepiped_width = symbols('parallelepiped_width')
parallelepiped_length = symbols('parallelepiped_length')
parallelepiped_height = symbols('parallelepiped_height')
geometric_buckling_squared = symbols('geometric_buckling_squared')

law = Eq(geometric_buckling_squared,
    (pi / parallelepiped_width)**2 + (pi / parallelepiped_length)**2 + (pi / parallelepiped_height)**2)

# This law is derived from geometric buckling definition (see geometric_buckling_from_neutron_flux.py),
# neutron flux laplacian in cartesian coordinates and boundary conditions.

# Unfortunately sympy does not support solving with complex boundary conditions so we simply check with known
# solution for the neutron flux:
# See [neutron flux for uniform parallelepiped](./neutron_flux_for_uniform_parallelepiped.py)
geometric_buckling_parallelepiped_squared = (parallelepiped_flux.width_constant**2 +
    parallelepiped_flux.length_constant**2 + parallelepiped_flux.height_constant**2)
geometric_buckling_parallelepiped_solved = geometric_buckling_parallelepiped_squared.subs({
    parallelepiped_flux.parallelepiped_width: parallelepiped_width,
    parallelepiped_flux.parallelepiped_length: parallelepiped_length,
    parallelepiped_flux.parallelepiped_height: parallelepiped_height
})
assert geometric_buckling_parallelepiped_solved == law.rhs

def print():
    return pretty(law, use_unicode=False)

@validate_input(parallelepiped_width_=units.length, parallelepiped_length_=units.length, parallelepiped_height_=units.length)
@validate_output(1 / units.length**2)
def calculate_geometric_buckling_squared(
    parallelepiped_width_: Quantity, parallelepiped_length_: Quantity, parallelepiped_height_: Quantity) -> Quantity:
    solved = solve(law, geometric_buckling_squared, dict=True)[0][geometric_buckling_squared]
    result_expr = solved.subs({
        parallelepiped_width: parallelepiped_width_,
        parallelepiped_length: parallelepiped_length_,
        parallelepiped_height: parallelepiped_height_})
    return expr_to_quantity(result_expr, 'parallelepiped_geometric_buckling_squared')
