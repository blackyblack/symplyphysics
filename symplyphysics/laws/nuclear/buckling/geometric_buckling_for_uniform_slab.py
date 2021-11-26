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

# Unfortunately sympy does not support solving with complex boundary conditions so we simply check with known
# solution for the neutron flux:
# See [neutron flux for uniform slab](./neutron_flux_for_uniform_slab.py)
geometric_buckling_slab_squared = slab_flux.axial_constant**2
geometric_buckling_slab_solved = geometric_buckling_slab_squared.subs(slab_flux.slab_width, slab_width)
assert geometric_buckling_slab_solved == law.rhs

def print():
    return pretty(law, use_unicode=False)

@validate_input(slab_width_=units.length)
@validate_output(1 / units.length**2)
def calculate_geometric_buckling_squared(slab_width_: Quantity) -> Quantity:
    solved = solve(law, geometric_buckling_squared, dict=True)[0][geometric_buckling_squared]
    result_expr = solved.subs(slab_width, slab_width_)
    return expr_to_quantity(result_expr, 'slab_geometric_buckling_squared')
