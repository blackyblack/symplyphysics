"""
Geometric buckling for uniform slab
===================================

Geometric buckling for the reactor in the shape of an infinite uniform slab of finite
thickness is a function of its thickness.

**Notes:**

#. Also see :ref:`Geometric buckling from neutron flux`.

**Links:**

#. `Wikipedia, derivable from third row of first table <https://en.wikipedia.org/wiki/Geometric_and_material_buckling#Geometric_Buckling>`__.
"""

from sympy import Eq, solve, pi
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.laws.nuclear.buckling import neutron_flux_for_uniform_slab as slab_flux

thickness = symbols.thickness
"""
:symbols:`thickness` of the slab.
"""

geometric_buckling = symbols.geometric_buckling
"""
:symbols:`geometric_buckling`.
"""

law = Eq(geometric_buckling, (pi / thickness)**2)
"""
:laws:symbol::

:laws:latex::
"""

# This law is derived from geometric buckling definition (see geometric_buckling_from_neutron_flux.py),
# neutron flux laplacian in cartesian coordinates and boundary condtitions.

# Unfortunately sympy does not support solving with complex boundary conditions so we simply check with known
# solution for the neutron flux:
# See [neutron flux for uniform slab](./neutron_flux_for_uniform_slab.py)
_geometric_buckling_slab_squared = slab_flux.axial_constant**2
_geometric_buckling_slab_solved = _geometric_buckling_slab_squared.subs(slab_flux.thickness,
    thickness)
assert _geometric_buckling_slab_solved == law.rhs

# TODO: derive from [parallelepiped law](./geometric_buckling_for_uniform_parallelepiped.py)

@validate_input(slab_width_=thickness)
@validate_output(geometric_buckling)
def calculate_geometric_buckling_squared(slab_width_: Quantity) -> Quantity:
    solved = solve(law, geometric_buckling, dict=True)[0][geometric_buckling]
    result_expr = solved.subs(thickness, slab_width_)
    return Quantity(result_expr)
