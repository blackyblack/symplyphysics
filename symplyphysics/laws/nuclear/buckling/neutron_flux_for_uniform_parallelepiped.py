"""
Neutron flux for uniform parallelepiped
=======================================

Neutron flux for a uniform rectangular parallelepiped reactor of side lengths
:math:`a, b, c` depends on the cartesian coordinates :math:`x, y, z`.
"""

from sympy import Eq, pi, cos
from sympy.vector import CoordSys3D
from symplyphysics import Quantity, units, symbols, clone_as_symbol
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_neutron_flux
from symplyphysics.laws.nuclear.buckling import neutron_flux_for_uniform_slab

dimension_factor = clone_as_symbol(symbols.neutron_flux, subscript="0")
"""
Dimension factor that appears as a coefficient in the solution to the :ref:`differential
equation <Diffusion equation from neutron flux>`. See :symbols:`neutron_flux`.
"""

x = clone_as_symbol(symbols.position, subscript="1")
"""
:symbols:`position` along the :math:`x`-axis.
"""

y = clone_as_symbol(symbols.position, subscript="2")
"""
:symbols:`position` along the :math:`y`-axis.
"""

z = clone_as_symbol(symbols.position, subscript="3")
"""
:symbols:`position` along the :math:`z`-axis.
"""

length = clone_as_symbol(symbols.length, subscript="1")
"""
:symbols:`length` along the :math:`x`-axis.
"""

width = clone_as_symbol(symbols.length, subscript="2")
"""
:symbols:`length` along the :math:`y`-axis.
"""

height = clone_as_symbol(symbols.length, subscript="3")
"""
:symbols:`length` along the :math:`z`-axis.
"""

neutron_flux = symbols.neutron_flux
"""
:symbols:`neutron_flux` at a point with coordinates :attr:`~x`, :attr:`~y`, :attr:`~z`.
"""

# These constants are being used for geometric buckling calculation
# See: [geometric buckling for uniform parallelepiped](geometric_buckling_for_uniform_parallelepiped.py)
_length_constant = pi / length
_width_constant = pi / width
_height_constant = pi / height

# derived the same way as uniform slab axial_constant
assert _length_constant == neutron_flux_for_uniform_slab._axial_constant.subs(  # pylint: disable=protected-access
    neutron_flux_for_uniform_slab.thickness, length)
assert _width_constant == neutron_flux_for_uniform_slab._axial_constant.subs(  # pylint: disable=protected-access
    neutron_flux_for_uniform_slab.thickness, width)
assert _height_constant == neutron_flux_for_uniform_slab._axial_constant.subs(  # pylint: disable=protected-access
    neutron_flux_for_uniform_slab.thickness, height)

law = Eq(
    neutron_flux,
    dimension_factor * cos(_width_constant * x) * cos(_length_constant * y) * cos(_height_constant * z))
"""
:laws:symbol::

:laws:latex::
"""

# Check the solution by passing the known neutron flux to the geometric_buckling_from_neutron_flux.
# Neutron flux is a function of x, y, z in the cartesian coordinates.

# Boundary conditions:
# - vacuum boundary condition: Ф(a + d) = Ф(ae) = Ф(b + d) = Ф(be) = Ф(c + d) = Ф(ce) = 0
# - finite flux condition: 0 <= Ф(x, y, z) < ∞
# - interface condition: the neutron flux and the normal component of the neutron current must be continuous
# - source condition: all neutrons flowing through the bounding area of the source must come from the neutron source

# define flux function in cylindrical coordinates as a function of cylinder radius and height
_cartesian_coordinates = CoordSys3D("_cartesian_coordinates")
# Make linter happy
_x = getattr(_cartesian_coordinates, "x")
_y = getattr(_cartesian_coordinates, "y")
_z = getattr(_cartesian_coordinates, "z")
_unit_length = Quantity(1, dimension=units.length)
neutron_flux_function_cartesian = law.subs({
    x: _x * _unit_length,
    y: _y * _unit_length,
    z: _z * _unit_length
})

_solved = geometric_buckling_from_neutron_flux.apply_neutron_flux_function(
    neutron_flux_function_cartesian.rhs)

# check with the derived law: Bg^2 = _width_constant**2 + _length_constant**2 + _height_constant**2
assert _solved.rhs == (_width_constant**2 + _length_constant**2 + _height_constant**2)


# There is no calculate() method. Neutron flux is usually being used internally to pass to other laws.
