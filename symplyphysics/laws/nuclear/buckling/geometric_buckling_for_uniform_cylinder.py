"""
Geometric buckling for uniform cylinder
=======================================

Geometric buckling of a uniform cylindrical reactor is a function of its radius
:math:`r` and height :math:`h`.

**Notes:**

#. Also see :ref:`Geometric buckling from neutron flux`.

**Links:**

#. `Wikipedia, second row in first table <https://en.wikipedia.org/wiki/Geometric_and_material_buckling#Geometric_Buckling>`__.
"""

from sympy import Eq, solve, pi
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.laws.nuclear.buckling import neutron_flux_for_uniform_cylinder as cylinder_flux

radius = symbols.radius
"""
:symbols:`radius` of the cylinder.
"""

height = symbols.height
"""
:symbols:`height` of the cylinder.
"""

geometric_buckling = symbols.geometric_buckling
"""
:symbols:`geometric_buckling`.
"""

law = Eq(geometric_buckling, (2.405 / radius)**2 + (pi / height)**2)
"""
:laws:symbol::

:laws:latex::
"""

# This law is derived from geometric buckling definition (see geometric_buckling_from_neutron_flux.py),
# neutron flux laplacian in cylindrical coordinates and boundary conditions.

# Unfortunately sympy does not support solving with complex boundary conditions so we simply check with known
# solution for the neutron flux:
# See [neutron flux for uniform cylinder](./neutron_flux_for_uniform_cylinder.py)
_geometric_buckling_cylinder_squared = cylinder_flux._radial_constant**2 + cylinder_flux._axial_constant**2
_geometric_buckling_cylinder_solved = _geometric_buckling_cylinder_squared.subs({
    cylinder_flux.radius: radius,
    cylinder_flux.height: height
})
assert (a := _geometric_buckling_cylinder_solved.evalf(7)) == (b := law.rhs.evalf(7)), (a, b)


@validate_input(cylinder_radius_=radius, cylinder_height_=height)
@validate_output(geometric_buckling)
def calculate_geometric_buckling_squared(cylinder_radius_: Quantity,
    cylinder_height_: Quantity) -> Quantity:
    solved = solve(law, geometric_buckling, dict=True)[0][geometric_buckling]
    result_expr = solved.subs({
        radius: cylinder_radius_,
        height: cylinder_height_
    })
    return Quantity(result_expr)
