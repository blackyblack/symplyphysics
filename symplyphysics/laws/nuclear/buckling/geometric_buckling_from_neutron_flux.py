from sympy import Function
from sympy.core.expr import Expr
from sympy.vector import Laplacian
from symplyphysics import (
    symbols, Eq, pretty, solve, units, simplify, expr_to_quantity
)
from symplyphysics.core.quantity_decorator import validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

# Description
## The quantity Bg^2 is called the geometrical buckling of the reactor and depends only on the geometry.
## This term is derived from the notion that the neutron flux distribution is somehow "buckled" in a homogeneous
## finite reactor. The neutron flux has more concave downward or "buckled" curvature (higher Bg^2) in a small
## reactor than in a large one.

## Law: Bg^2 = -1 * Δ(Ф(x)) / Ф(x)
## Where:
## Δ - Laplacian operator.
## Ф(x) (neutron flux density) - number of neutrons crossing through some arbitrary cross-sectional unit area in all
##   directions per unit time.
## x - coordinates for the neutron flux density.
## Bg^2 - geometric buckling.

#TODO: try to add proper dimensions
flux_position = symbols("flux_position")
neutron_flux = symbols("neutron_flux", cls=Function)
geometric_buckling_squared = Symbol("geometric_buckling_squared", 1 / units.length**2)

# neutron_flux_function should be a function on CoordSys3D, eg:
#   spherical_coordinates = CoordSys3D("spherical_coordinates", transformation="spherical")
#   neutron_flux_function(spherical_coordinates.r)
law = Eq(geometric_buckling_squared,
    -1 * Laplacian(neutron_flux(flux_position)) / neutron_flux(flux_position))

def print(expr: Expr) -> str:
    symbols = [geometric_buckling_squared]
    return pretty(to_printable(expr, symbols), use_unicode=False)

# This is mostly internal method for solving derived laws. End-user is unlikely to define
# neutron flux function and pass it here. See geometric_buckling_for_uniform_sphere.py as an example
# of such derived law.
# neutron_flux_function_ should be a function on CoordSys3D
def apply_neutron_flux_function(neutron_flux_function_: Function) -> Expr:
    applied_law = law.subs(neutron_flux(flux_position), neutron_flux_function_)
    return simplify(applied_law.doit())

# neutron_flux_function_ should be a function on CoordSys3D
# neutron_flux_function_ geometry should be defined with Quantity, eg width.dimension == units.length
@validate_output_symbol(geometric_buckling_squared)
def calculate_geometric_buckling_squared(neutron_flux_function_: Function) -> Quantity:
    result_expr = apply_neutron_flux_function(neutron_flux_function_)
    result_buckling_expr = solve(result_expr, geometric_buckling_squared, dict=True)[0][geometric_buckling_squared]
    return expr_to_quantity(result_buckling_expr)
