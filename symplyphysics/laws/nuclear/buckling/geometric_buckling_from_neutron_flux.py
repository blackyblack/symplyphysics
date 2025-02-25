"""
Geometric buckling from neutron flux
====================================

Geometric buckling can be found if the neutron flux density is known as a function of
position. The smaller the reactor is, more "buckled" the curvature of the neutron flux
is and the higher geometric buckling of the reactor is.

**Links:**

#. `Wikipedia, see fourth equation <https://en.wikipedia.org/wiki/Geometric_and_material_buckling#Derivation>`__.
"""

from sympy import (Eq, solve, Expr, simplify, Equality)
from sympy.vector import Laplacian
from symplyphysics import (SI, units, Quantity, validate_output, symbols, clone_as_function,
    Function)
from symplyphysics.core.dimensions import collect_quantity_factor_and_dimension
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.nuclear import diffusion_equation_from_neutron_flux as diffusion_equation

position = symbols.position
"""
:symbols:`position`.
"""

neutron_flux = clone_as_function(symbols.neutron_flux, [position])
"""
:symbols:`neutron_flux` as a function of :attr:`~position`
"""

geometric_buckling = symbols.geometric_buckling
"""
:symbols:`geometric_buckling`.
"""

neutron_flux_laplacian = Function("Laplace(Phi)", [position],
    1 / (units.length**4 * units.time),
    display_latex="\\nabla^{2} \\Phi")
"""
Laplacian of the :attr:`~neutron_flux` as a function of :attr:`~position`.
"""

# As Laplacian is a second derivative over space coordinates (x, y, z), resulting dimension should be
# original dimension / units.length**2
# pylint: disable-next=comparison-with-callable
assert neutron_flux_laplacian.dimension == diffusion_equation.neutron_flux_laplacian.dimension

neutron_flux_laplacian_definition = Eq(neutron_flux_laplacian(position),
    Laplacian(neutron_flux(position)),
    evaluate=False)

# neutron_flux_function should be a function on CoordSys3D, eg:
#   spherical_coordinates = CoordSys3D("spherical_coordinates", transformation="spherical")
#   neutron_flux_function(spherical_coordinates.r)
law = Eq(geometric_buckling, -1 * neutron_flux_laplacian(position) / neutron_flux(position))
"""
:laws:symbol::

:laws:latex::
"""

# Check laplacian definition is the same as in diffusion equation

_diffusion_equation_laplacian = diffusion_equation.neutron_flux_laplacian_definition.rhs.subs(
    diffusion_equation.neutron_flux(diffusion_equation.position), neutron_flux(position))
assert expr_equals(_diffusion_equation_laplacian, neutron_flux_laplacian_definition.rhs)


# neutron_flux_function_ should be a function on CoordSys3D
# This is an exact copy from 'diffusion_equation_from_neutron_flux'
def apply_neutron_flux_function(neutron_flux_function_: Expr) -> Equality:
    # Manually divide to unit_length to get Laplacian dimension. CoordSys3D coordinates are dimensionless, hence
    # Laplacian cannot properly calculate resulting dimension.
    unit_length = Quantity(1, dimension=units.length)
    neutron_flux_laplacian_eval = neutron_flux_laplacian_definition.rhs.subs(
        neutron_flux(position), neutron_flux_function_).doit() / unit_length**2
    applied_law = law.subs(neutron_flux_laplacian(position), neutron_flux_laplacian_eval)
    applied_law = applied_law.subs(neutron_flux(position), neutron_flux_function_)
    return simplify(applied_law)


# neutron_flux_function_ should be a function on CoordSys3D
# neutron_flux_function_ geometry should be defined with Quantity, eg width.dimension == units.length
@validate_output(geometric_buckling)
def calculate_geometric_buckling_squared(neutron_flux_function_: Expr, radius: Expr) -> Quantity:
    (_, dimension) = collect_quantity_factor_and_dimension(neutron_flux_function_.subs(radius, 1))
    assert SI.get_dimension_system().equivalent_dims(dimension, neutron_flux.dimension)

    result_expr = apply_neutron_flux_function(neutron_flux_function_)
    result_buckling_expr = solve(result_expr, geometric_buckling, dict=True)[0][geometric_buckling]
    return Quantity(result_buckling_expr)
