"""
Neutron flux for uniform slab
=============================

Neutron flux for a reactor in the shape of an infinite slab of finite thickness depends
on the orthogonal distance to the plane of symmetry inside the slab.

**Links:**

#. `NuclearPower, see end of page <https://www.nuclear-power.com/nuclear-power/reactor-physics/neutron-diffusion-theory/infinite-slab-reactor/>`__.
"""

from sympy import Eq, pi, cos
from sympy.vector import CoordSys3D
from symplyphysics import Quantity, units, symbols, clone_as_function, clone_as_symbol
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_neutron_flux

dimension_factor = clone_as_symbol(symbols.neutron_flux, subscript="0")
"""
Dimension factor that appears as a coefficient in the solution to the :ref:`differential
equation <Diffusion equation from neutron flux>`. See :symbols:`neutron_flux`.
"""

distance = symbols.orthogonal_distance
"""
:symbols:`orthogonal_distance` to the central plane of the slab.
"""

thickness = symbols.thickness
"""
:symbols:`thickness` of the slab.
"""

neutron_flux = symbols.neutron_flux
"""
:symbols:`neutron_flux` at a :attr:`~distance` from the central plane of the slab.
"""

# This constant is being used for geometric buckling calculation
# See: [geometric buckling for uniform slab](geometric_buckling_for_uniform_slab.py)
_axial_constant = pi / thickness

law = Eq(neutron_flux, dimension_factor * cos(_axial_constant * distance))
"""
:laws:symbol::

:laws:latex::
"""

# Check the solution by passing the known neutron flux to the geometric_buckling_from_neutron_flux.
# Neutron flux is a function of x coordinate in the cartesian coordinates.

# Boundary conditions:
# - vacuum boundary condition: Ф(a / 2 + d) = Ф(ae / 2) = 0
# - finite flux condition: 0 <= Ф(x) < ∞
# - interface condition: the neutron flux and the normal component of the neutron current must be continuous
# - source condition: all neutrons flowing through the bounding area of the source must come from the neutron source

# define flux function in cartesian coordinates as a function of x coordinate
_cartesian_coordinates = CoordSys3D("_cartesian_coordinates")
# Make linter happy
_x = getattr(_cartesian_coordinates, "x")
_unit_length = Quantity(1, dimension=units.length)
_neutron_flux_function_cartesian = law.subs(distance, _x * _unit_length)

_solved = geometric_buckling_from_neutron_flux.apply_neutron_flux_function(
    _neutron_flux_function_cartesian.rhs)

# check with the derived law: Bg^2 = _axial_constant**2
assert _solved.rhs == _axial_constant**2


# There is no calculate() method. Neutron flux is usually being used internally to pass to other laws.
