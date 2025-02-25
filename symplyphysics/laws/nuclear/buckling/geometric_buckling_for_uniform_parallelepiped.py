"""
Geometric buckling for uniform parallelepiped
=============================================

Geometric buckling of a uniform parallelepiped reactor is a function of its side lengths
:math:`a, b, c`.

**Notes:**

#. Also see :ref:`Geometric buckling from neutron flux`.

**Links:**

#. `Wikipedia, third row in first table <https://en.wikipedia.org/wiki/Geometric_and_material_buckling#Geometric_Buckling>`__.
"""

from sympy import (Eq, solve, pi)
from symplyphysics import (Quantity, validate_input, validate_output, symbols, clone_as_symbol)
from symplyphysics.laws.nuclear.buckling import neutron_flux_for_uniform_parallelepiped as parallelepiped_flux

length = clone_as_symbol(symbols.length, subscript="1")
"""
:symbols:`length` of the first side of the parallelepiped.
"""

width = clone_as_symbol(symbols.length, subscript="2")
"""
:symbols:`length` of the second side of the parallelepiped.
"""

height = clone_as_symbol(symbols.length, subscript="3")
"""
:symbols:`length` of the third side of the parallelepiped.
"""

geometric_buckling = symbols.geometric_buckling
"""
:symbols:`geometric_buckling`.
"""

law = Eq(geometric_buckling, (pi / length)**2 + (pi / width)**2 + (pi / height)**2)
"""
:laws:symbol::

:laws:latex::
"""

# This law is derived from geometric buckling definition (see geometric_buckling_from_neutron_flux.py),
# neutron flux laplacian in cartesian coordinates and boundary conditions.

# Unfortunately sympy does not support solving with complex boundary conditions so we simply check with known
# solution for the neutron flux:
# See [neutron flux for uniform parallelepiped](./neutron_flux_for_uniform_parallelepiped.py)
_geometric_buckling_parallelepiped_squared = (parallelepiped_flux.width_constant**2 +
    parallelepiped_flux.length_constant**2 + parallelepiped_flux.height_constant**2)
_geometric_buckling_parallelepiped_solved = _geometric_buckling_parallelepiped_squared.subs({
    parallelepiped_flux.width: width,
    parallelepiped_flux.length: length,
    parallelepiped_flux.height: height
})
assert _geometric_buckling_parallelepiped_solved == law.rhs


@validate_input(parallelepiped_width_=width,
    parallelepiped_length_=length,
    parallelepiped_height_=height)
@validate_output(geometric_buckling)
def calculate_geometric_buckling_squared(parallelepiped_width_: Quantity,
    parallelepiped_length_: Quantity, parallelepiped_height_: Quantity) -> Quantity:
    solved = solve(law, geometric_buckling, dict=True)[0][geometric_buckling]
    result_expr = solved.subs({
        width: parallelepiped_width_,
        length: parallelepiped_length_,
        height: parallelepiped_height_
    })
    return Quantity(result_expr)
