from sympy import pi
from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity, filter_map_zeroes
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

# Boundary conditions:
# - vacuum boundary condition: Ф(R + d) =  Ф(Re) = 0
# - finite flux condition: 0 <= Ф(r) < ∞
# - interface condition: the neutron flux and the normal component of the neutron current must be continuous
# - source condition: all neutrons flowing through the bounding area of the source must come from the neutron source
# - albedo boundary condition: Ф(Ralbedo) = 0

# Unfortunately sympy does not support solving with complex boundary conditions so we simply check with known
# solution for the neutron flux:
# See [neutron flux for uniform sphere](./neutron_flux_for_uniform_sphere.py)
vacuum_boundary_condition = Eq(sphere_flux.neutron_flux_function(sphere_radius), 0)
vacuum_boundary_neutron_flux_function = sphere_flux.law.subs(sphere_flux.distance_from_center, sphere_radius)

geometric_buckling_solved = solve([vacuum_boundary_neutron_flux_function, vacuum_boundary_condition],
    (sphere_flux.neutron_flux_function(sphere_radius), sphere_flux.geometric_buckling_sphere), dict=True)

# first solution is zero. We do not need it due to boundary conditions.
geometric_buckling_solved = filter_map_zeroes(sphere_flux.geometric_buckling_sphere, geometric_buckling_solved)

geometric_buckling_solution = geometric_buckling_solved[0][sphere_flux.geometric_buckling_sphere]
assert geometric_buckling_solution**2 == law.rhs

def print():
    return pretty(law, use_unicode=False)

@validate_input(sphere_radius_=units.length)
@validate_output(1 / units.length**2)
def calculate_geometric_buckling_squared(sphere_radius_: Quantity) -> Quantity:
    solved = solve(law, geometric_buckling_squared, dict=True)[0][geometric_buckling_squared]
    result_expr = solved.subs(sphere_radius, sphere_radius_)
    return expr_to_quantity(result_expr, 'sphere_geometric_buckling_squared')
