"""
Geometric buckling for uniform sphere
=====================================

Geometric buckling for a uniform spherical reactor is a function of its radius.

**Links:**

#. `Wikipedia, first row in first table <https://en.wikipedia.org/wiki/Geometric_and_material_buckling#Geometric_Buckling>`__.
"""

from sympy import Eq, solve, pi
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.laws.nuclear.buckling import neutron_flux_for_uniform_sphere as sphere_flux

radius = symbols.radius
"""
:symbols:`radius` of the sphere.
"""

geometric_buckling = symbols.geometric_buckling
"""
:symbols:`geometric_buckling`.
"""

law = Eq(geometric_buckling, (pi / radius)**2)
"""
:laws:symbol::

:laws:latex::
"""

# This law is derived from geometric buckling definition (see geometric_buckling_from_neutron_flux.py),
# neutron flux laplacian in spherical coordinates and boundary condtitions.

# Unfortunately sympy does not support solving with complex boundary conditions so we simply check with known
# solution for the neutron flux:
# See [neutron flux for uniform sphere](./neutron_flux_for_uniform_sphere.py)
_geometric_buckling_sphere_squared = sphere_flux.radial_constant**2
_geometric_buckling_sphere_flux_solved = _geometric_buckling_sphere_squared.subs(
    sphere_flux.radius, radius)
assert _geometric_buckling_sphere_flux_solved == law.rhs


@validate_input(sphere_radius_=radius)
@validate_output(geometric_buckling)
def calculate_geometric_buckling_squared(sphere_radius_: Quantity) -> Quantity:
    solved = solve(law, geometric_buckling, dict=True)[0][geometric_buckling]
    result_expr = solved.subs(radius, sphere_radius_)
    return Quantity(result_expr)
