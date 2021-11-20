from sympy import pi, cos
from sympy.vector import CoordSys3D
from symplyphysics import (
    symbols, Function, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity, filter_map_zeroes
)
from symplyphysics.laws.nuclear import geometric_buckling_from_neutron_flux

# Description
## Geometric buckling for the reactor in the shape of a slab of physical width a in the x-direction
## and infinite in the y- and z-directions.

## Law: Bg^2slab = (Pi / a)^2
## Where:
## Pi - Pi constant.
## a - slab width.
## Bg^2slab - geometric buckling for slab.
##   See [geometric buckling](./geometric_buckling_from_neutron_flux.py) implementation.

slab_width = symbols('slab_width')
geometric_buckling = symbols('geometric_buckling')

law = Eq(geometric_buckling, (pi / slab_width)**2)

# Derive the same law from diffusion equation and laplacian in rectangle coordinates

# this solution should be derived from the boundary conditions:
# - vacuum boundary condition: Ф(a / 2 + d) =  Ф(ae / 2) = 0
# - finite flux condition: 0 <= Ф(x) < ∞
# - interface condition: the neutron flux and the normal component of the neutron current must be continuous
# - source condition: all neutrons flowing through the bounding area of the source must come from the neutron source
# unfortunately sympy does not support solving with complex boundary conditions so we simply check with one of
# the known solutions:
neutron_flux_power_constant = symbols('C1', constant=True)
# we know that it should be constant but we pass it as a function to the dsolve()
geometric_buckling_sqrt = symbols('geometric_buckling_sqrt')
geometric_buckling_function = symbols('geometric_buckling_function', cls = Function)
distance_from_center = symbols('distance_from_center')
neutron_flux_function = symbols('neutron_flux_function', cls = Function)

# this is an expected solution for the slab reactor
neutron_flux_function_solution = Eq(neutron_flux_function(distance_from_center),
    neutron_flux_power_constant * cos(geometric_buckling_sqrt * distance_from_center))

vacuum_boundary_condition = Eq(neutron_flux_function(slab_width / 2), 0)
vacuum_boundary_neutron_flux_function = neutron_flux_function_solution.subs(distance_from_center, slab_width / 2)

geometric_buckling_sqrt_solved = solve([vacuum_boundary_neutron_flux_function, vacuum_boundary_condition],
    (neutron_flux_function(slab_width / 2), geometric_buckling_sqrt), dict=True)

# we only accept first solution because other solutions violate finite flux condition
geometric_buckling_sqrt_solution = geometric_buckling_sqrt_solved[0][geometric_buckling_sqrt]

derived_geometric_buckling = geometric_buckling_sqrt_solution**2
assert derived_geometric_buckling == law.rhs

# check again by passing calculated geometric buckling value as a parameter for the neutron flux
# that is used to calculate geometric buckling squared.
# This check is redundant - better move it to more general law checks.

# define flux function in rectangular coordinates as a function of x position
rectangular_coordinates = CoordSys3D('rectangular_coordinates')
neutron_flux_function_rectangular = neutron_flux_function_solution.subs({
    distance_from_center: rectangular_coordinates.x,
    geometric_buckling_sqrt: geometric_buckling_sqrt_solution})

solved = geometric_buckling_from_neutron_flux.calculate_geometric_buckling(neutron_flux_function_rectangular.rhs)

# check with the derived law
assert solved.rhs == law.rhs

def print():
    return pretty(law, use_unicode=False)

@validate_input(slab_width_=units.length)
@validate_output(1 / units.length**2)
def calculate_geometric_buckling_squared(slab_width_: Quantity) -> Quantity:
    solved = solve(law, geometric_buckling, dict=True)[0][geometric_buckling]
    result_expr = solved.subs(slab_width, slab_width_)
    return expr_to_quantity(result_expr, 'slab_geometric_buckling_squared')
