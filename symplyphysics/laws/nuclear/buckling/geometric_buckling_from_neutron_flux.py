from sympy import (Eq, solve, Expr, symbols, simplify)
from sympy.vector import Laplacian
from symplyphysics import (SI, Function, units, expr_to_quantity, Quantity, Symbol,
    print_expression, validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.symbols.quantities import collect_factor_and_dimension
from symplyphysics.laws.nuclear import diffusion_equation_from_neutron_flux as diffusion_equation

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

# Position is a free variable of a function - do not specify its dimension
flux_position = symbols("flux_position")
neutron_flux = Function("neutron_flux", 1 / units.length**2 / units.time)
geometric_buckling_squared = Symbol("geometric_buckling_squared", 1 / units.length**2)
neutron_flux_laplacian = Function("neutron_flux_laplacian", 1 / units.length**4 / units.time)

# As Laplacian is a second derivative over space coordinates (x, y, z), resulting dimension should be
# original dimension / units.length**2
assert neutron_flux_laplacian.dimension == diffusion_equation.neutron_flux_laplacian.dimension

neutron_flux_laplacian_definition = Eq(neutron_flux_laplacian(flux_position),
    Laplacian(neutron_flux(flux_position)),
    evaluate=False)

# neutron_flux_function should be a function on CoordSys3D, eg:
#   spherical_coordinates = CoordSys3D("spherical_coordinates", transformation="spherical")
#   neutron_flux_function(spherical_coordinates.r)
law = Eq(geometric_buckling_squared,
    -1 * neutron_flux_laplacian(flux_position) / neutron_flux(flux_position))

# Check laplacian definition is the same as in diffusion equation

diffusion_equation_laplacian = diffusion_equation.neutron_flux_laplacian_definition.rhs.subs(
    diffusion_equation.neutron_flux(diffusion_equation.flux_position), neutron_flux(flux_position))
assert expr_equals(diffusion_equation_laplacian, neutron_flux_laplacian_definition.rhs)


def print_law() -> str:
    return print_expression(law)


# neutron_flux_function_ should be a function on CoordSys3D
# This is an exact copy from 'diffusion_equation_from_neutron_flux'
def apply_neutron_flux_function(neutron_flux_function_: Expr) -> Expr:
    # Manually divide to unit_length to get Laplacian dimension. CoordSys3D coordinates are dimensionless, hence
    # Laplacian cannot properly calculate resulting dimension.
    unit_length = Quantity(1, dimension=units.length)
    neutron_flux_laplacian_eval = neutron_flux_laplacian_definition.rhs.subs(
        neutron_flux(flux_position), neutron_flux_function_).doit() / unit_length**2
    applied_law = law.subs(neutron_flux_laplacian(flux_position), neutron_flux_laplacian_eval)
    applied_law = applied_law.subs(neutron_flux(flux_position), neutron_flux_function_)
    return simplify(applied_law)


# neutron_flux_function_ should be a function on CoordSys3D
# neutron_flux_function_ geometry should be defined with Quantity, eg width.dimension == units.length
@validate_output(geometric_buckling_squared)
def calculate_geometric_buckling_squared(neutron_flux_function_: Expr) -> Quantity:
    (_, neutron_flux_dimension) = collect_factor_and_dimension(neutron_flux_function_)
    assert SI.get_dimension_system().equivalent_dims(neutron_flux_dimension, neutron_flux.dimension)

    result_expr = apply_neutron_flux_function(neutron_flux_function_)
    result_buckling_expr = solve(result_expr, geometric_buckling_squared,
        dict=True)[0][geometric_buckling_squared]
    return expr_to_quantity(result_buckling_expr)
