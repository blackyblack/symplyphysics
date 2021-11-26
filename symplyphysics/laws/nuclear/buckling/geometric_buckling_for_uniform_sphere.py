from sympy import pi
from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)
from symplyphysics.laws.nuclear.buckling import neutron_flux_for_uniform_sphere as sphere_flux

# Description
## Geometric buckling for the uniform spherical reactor. The spherical reactor is situated in spherical
## geometry at the origin of coordinates.

## Law: Bg^2sphere = (Pi / R)^2
## Where:
## Pi - Pi constant.
## R - sphere radius.
## Bg^2sphere - squared geometric buckling for sphere.
##   See [geometric buckling](./geometric_buckling_from_neutron_flux.py) implementation.

sphere_radius = symbols('sphere_radius')
geometric_buckling_squared = symbols('geometric_buckling_squared')

law = Eq(geometric_buckling_squared, (pi / sphere_radius)**2)

# This law is derived from geometric buckling definition (see geometric_buckling_from_neutron_flux.py),
# neutron flux laplacian in spherical coordinates and boundary condtitions.

# Unfortunately sympy does not support solving with complex boundary conditions so we simply check with known
# solution for the neutron flux:
# See [neutron flux for uniform sphere](./neutron_flux_for_uniform_sphere.py)
geometric_buckling_sphere_squared = sphere_flux.radial_constant**2
geometric_buckling_sphere_flux_solved = geometric_buckling_sphere_squared.subs(sphere_flux.sphere_radius, sphere_radius)
assert geometric_buckling_sphere_flux_solved == law.rhs

def print():
    return pretty(law, use_unicode=False)

@validate_input(sphere_radius_=units.length)
@validate_output(1 / units.length**2)
def calculate_geometric_buckling_squared(sphere_radius_: Quantity) -> Quantity:
    solved = solve(law, geometric_buckling_squared, dict=True)[0][geometric_buckling_squared]
    result_expr = solved.subs(sphere_radius, sphere_radius_)
    return expr_to_quantity(result_expr, 'sphere_geometric_buckling_squared')
