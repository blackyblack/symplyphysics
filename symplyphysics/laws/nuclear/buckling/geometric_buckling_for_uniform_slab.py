from sympy import (Eq, solve, pi)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input, validate_output)
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

slab_width = Symbol("slab_width", units.length)
geometric_buckling_squared = Symbol("geometric_buckling_squared", 1 / units.area)

law = Eq(geometric_buckling_squared, (pi / slab_width)**2)

# This law is derived from geometric buckling definition (see geometric_buckling_from_neutron_flux.py),
# neutron flux laplacian in cartesian coordinates and boundary condtitions.

# Unfortunately sympy does not support solving with complex boundary conditions so we simply check with known
# solution for the neutron flux:
# See [neutron flux for uniform slab](./neutron_flux_for_uniform_slab.py)
geometric_buckling_slab_squared = slab_flux.axial_constant**2
geometric_buckling_slab_solved = geometric_buckling_slab_squared.subs(slab_flux.slab_width,
    slab_width)
assert geometric_buckling_slab_solved == law.rhs


def print_law() -> str:
    return print_expression(law)


@validate_input(slab_width_=slab_width)
@validate_output(geometric_buckling_squared)
def calculate_geometric_buckling_squared(slab_width_: Quantity) -> Quantity:
    solved = solve(law, geometric_buckling_squared, dict=True)[0][geometric_buckling_squared]
    result_expr = solved.subs(slab_width, slab_width_)
    return expr_to_quantity(result_expr)
