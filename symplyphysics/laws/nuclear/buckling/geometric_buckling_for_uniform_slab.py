from sympy import pi
from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)
from symplyphysics.laws.nuclear.buckling import neutron_flux_for_uniform_slab as slab_flux

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
geometric_buckling_squared = symbols('geometric_buckling_squared')

law = Eq(geometric_buckling_squared, (pi / slab_width)**2)

# This law is derived from geometric buckling definition (see geometric_buckling_from_neutron_flux.py),
# neutron flux laplacian in cartesian coordinates and boundary condtitions.

# Boundary conditions:
# - vacuum boundary condition: Ф(a / 2 + d) =  Ф(ae / 2) = 0
# - finite flux condition: 0 <= Ф(x) < ∞
# - interface condition: the neutron flux and the normal component of the neutron current must be continuous
# - source condition: all neutrons flowing through the bounding area of the source must come from the neutron source

# Unfortunately sympy does not support solving with complex boundary conditions so we simply check with known
# solution for the neutron flux:
# See [neutron flux for uniform slab](./neutron_flux_for_uniform_slab.py)
vacuum_boundary_condition = Eq(slab_flux.neutron_flux_function(slab_width / 2), 0)
vacuum_boundary_neutron_flux_function = slab_flux.law.subs(slab_flux.distance_from_center, slab_width / 2)

geometric_buckling_solved = solve([vacuum_boundary_neutron_flux_function, vacuum_boundary_condition],
    (slab_flux.neutron_flux_function(slab_width / 2), slab_flux.geometric_buckling_slab), dict=True)

# we only accept first solution because other solutions violate finite flux condition
geometric_buckling_solution = geometric_buckling_solved[0][slab_flux.geometric_buckling_slab]
assert geometric_buckling_solution**2 == law.rhs

def print():
    return pretty(law, use_unicode=False)

@validate_input(slab_width_=units.length)
@validate_output(1 / units.length**2)
def calculate_geometric_buckling_squared(slab_width_: Quantity) -> Quantity:
    solved = solve(law, geometric_buckling_squared, dict=True)[0][geometric_buckling_squared]
    result_expr = solved.subs(slab_width, slab_width_)
    return expr_to_quantity(result_expr, 'slab_geometric_buckling_squared')
