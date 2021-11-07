from sympy import pi, sin, simplify
from sympy.vector import CoordSys3D, laplacian
from symplyphysics import (
    symbols, Function, Eq, pretty, solve, dsolve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Geometric buckling for the uniform spherical reactor. The spherical reactor is situated in spherical
## geometry at the origin of coordinates.

## Law: Bg^2sphere = (Pi / R)^2
## Where:
## Pi - Pi constant.
## R - sphere radius.
## Bg^2sphere - geometric buckling.
##   See [geometric buckling](./geometric_buckling_from_neutron_flux.py)

sphere_radius = symbols('sphere_radius')
geometric_buckling = symbols('geometric_buckling')

law = Eq(geometric_buckling, (pi / sphere_radius)**2)

# Derive the same law from diffusion equation and laplacian in spherical coordinates

distance_from_center = symbols('distance_from_center')
neutron_flux_function = symbols('neutron_flux_function', cls = Function)

# derive from spherical coordinates Laplacian
spherical_coordinates = CoordSys3D('spherical_coordinates', transformation='spherical')
# flux does not depend on angles so it is a funcion of sphere radius only
flux_spherical_laplacian_with_theta_angle = laplacian(neutron_flux_function(spherical_coordinates.r))
# laplacian still contains an angle but we do not need it - flux should not depend on angle rotation
flux_spherical_laplacian = flux_spherical_laplacian_with_theta_angle.subs(sin(spherical_coordinates.theta), 1)
# use distance_from_center as parameter to allow sympy to process this expression
flux_laplacian_function = flux_spherical_laplacian.subs(spherical_coordinates.r, distance_from_center)

# see geometric buckling definition
geometric_buckling_from_diffusion_law = simplify(Eq(geometric_buckling,
    -1 * flux_laplacian_function / neutron_flux_function(distance_from_center)))

# this solution should be derived from the boundary conditions:
# - vacuum boundary condition: Ф(R + d) =  Ф(Re) = 0
# - finite flux condition: 0 <= Ф(r) < ∞
# - interface condition: the neutron flux and the normal component of the neutron current must be continuous
# - source condition: all neutrons flowing through the bounding area of the source must come from the neutron source
# - albedo boundary condition: Ф(Ralbedo) = 0
# unfortunately sympy does not support solving with complex boundary conditions so we simply check with one of
# the known solutions:
neutron_flux_power_constant = symbols('C1', constant=True)
# we know that it should be constant but we pass it as a function to the dsolve()
geometric_buckling_sqrt = symbols('geometric_buckling_sqrt')
geometric_buckling_function = symbols('geometric_buckling_function', cls = Function)

# this is an expected solution for the spherical reactor
neutron_flux_function_solution = Eq(neutron_flux_function(distance_from_center),
    neutron_flux_power_constant * sin(geometric_buckling_sqrt * distance_from_center) / distance_from_center)

vacuum_boundary_condition = Eq(neutron_flux_function(sphere_radius), 0)
vacuum_boundary_neutron_flux_function = neutron_flux_function_solution.subs(distance_from_center, sphere_radius)

geometric_buckling_sqrt_solved = solve([vacuum_boundary_neutron_flux_function, vacuum_boundary_condition],
    (neutron_flux_function(sphere_radius), geometric_buckling_sqrt), dict=True)

# first solution is zero. We do not need it due to boundary conditions.
# TODO: add filters for zero and negative solutions
assert geometric_buckling_sqrt_solved[0][geometric_buckling_sqrt] == 0

geometric_buckling_sqrt_solution = geometric_buckling_sqrt_solved[1][geometric_buckling_sqrt]
derived_geometric_buckling = geometric_buckling_sqrt_solution**2
assert derived_geometric_buckling == law.rhs

# check again by passing calculated geometric buckling value as a parameter for the neutron flux
# that is used to calculate geometric buckling squared.
# This check is redundant - better move it to more general law checks.

neutron_flux_function_spherical = neutron_flux_function_solution.subs(
    geometric_buckling_sqrt, geometric_buckling_sqrt_solution)
applied_derived_law = geometric_buckling_from_diffusion_law.subs({
    neutron_flux_function(distance_from_center): neutron_flux_function_spherical.rhs,
    geometric_buckling: geometric_buckling_function(distance_from_center)})
solved = dsolve(applied_derived_law, geometric_buckling_function(distance_from_center))

# check with the derived law
assert solved.rhs == law.rhs

def print():
    return pretty(law, use_unicode=False)

@validate_input(sphere_radius_=units.length)
@validate_output(1 / units.length**2)
def calculate_geometric_buckling_squared(sphere_radius_: Quantity) -> Quantity:
    solved = solve(law, geometric_buckling, dict=True)[0][geometric_buckling]
    result_expr = solved.subs({
        sphere_radius: sphere_radius_})
    return expr_to_quantity(result_expr, 'sphere_geometric_buckling_squared')
