"""
Neutron flux for uniform cylinder
=================================

Neutron flux for the uniform cylindrical reactor of radius :math:`r_0` and height
:math:`h_0` depends on the radial distance :math:`r` and axial coordinate :math:`z`. A
cylindrical coordinate system is assumed, such that the cylinder's axis of rotation is
the :math:`z`-axis.

**Links:**

#. `NuclearPower, see end of page <https://www.nuclear-power.com/nuclear-power/reactor-physics/neutron-diffusion-theory/finite-cylindrical-reactor/>`__.
"""

from sympy import Eq, pi, cos
from sympy.vector import CoordSys3D
from sympy.functions.special.bessel import besselj
from symplyphysics import Quantity, units, symbols, clone_as_symbol
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_neutron_flux
from symplyphysics.laws.nuclear.buckling import neutron_flux_for_uniform_slab

dimension_factor = clone_as_symbol(symbols.neutron_flux, subscript="0")
"""
Dimension factor that appears as a coefficient in the solution to the :ref:`differential
equation <Diffusion equation from neutron flux>`. See :symbols:`neutron_flux`.
"""

radial_distance = symbols.distance_to_axis
"""
:symbols:`distance_to_axis` within the cylindrical coordinate system of the cylinder.
"""

axial_coordinate = symbols.height
"""
:symbols:`height` within the cylindrical coordinate system of the cylinder.
"""

radius = clone_as_symbol(symbols.radius, subscript="0")
"""
:symbols:`radius` of the cylinder.
"""

height = clone_as_symbol(symbols.height, subscript="0")
"""
:symbols:`height` of the cylinder.
"""

neutron_flux = symbols.neutron_flux
"""
:symbols:`neutron_flux` at a point with coordinates :attr:`~radial_distance` and
:attr:`~axial_coordinate`.
"""

# These constants are being used for geometric buckling calculation
# See: [geometric buckling for uniform cylinder](geometric_buckling_for_uniform_cylinder.py)
_radial_constant = 2.405 / radius
_axial_constant = pi / height

# derived the same way as uniform slab _axial_constant
assert _axial_constant == neutron_flux_for_uniform_slab._axial_constant.subs(  # pylint: disable=protected-access
    neutron_flux_for_uniform_slab.thickness, height)

law = Eq(
    neutron_flux,
    dimension_factor * besselj(0, _radial_constant * radial_distance) *
    cos(_axial_constant * axial_coordinate))
"""
:laws:symbol::

:laws:latex::
"""

# Check the solution by passing the known neutron flux to the geometric_buckling_from_neutron_flux.
# Neutron flux is a function of radius and height in the cylindrical coordinates.

# Boundary conditions:
# - vacuum boundary condition: Ф(R + d) = Ф(Re) = Ф(a + d) = Ф(ae) = 0
# - finite flux condition: 0 <= Ф(r, x) < ∞
# - interface condition: the neutron flux and the normal component of the neutron current must be continuous
# - source condition: all neutrons flowing through the bounding area of the source must come from the neutron source
# - albedo boundary condition: Ф(Ralbedo) = 0

# _radial_constant is the solution of the Bessel function J0, with a condition the neutron flux cannot
# have negative values (finite flux condition) and with zero flux boundary condition

# Define flux function in cylindrical coordinates as a function of cylinder radius and height

# CoordinateSystem class does not work here, because Laplacian obtains coordinate system from
# the provided scalar field (neutron_flux function)
_cylindrical_coordinates = CoordSys3D("_cylindrical_coordinates", transformation="cylindrical")
# Make linter happy
_r = getattr(_cylindrical_coordinates, "r")
_z = getattr(_cylindrical_coordinates, "z")
_unit_length = Quantity(1, dimension=units.length)
_neutron_flux_function_cylindrical = law.subs({
    radial_distance: _r * _unit_length,
    axial_coordinate: _z * _unit_length
})

_solved = geometric_buckling_from_neutron_flux.apply_neutron_flux_function(
    _neutron_flux_function_cylindrical.rhs)

# check with the derived law: Bg^2 = _radial_constant**2 + _axial_constant**2
# limit decimals to bypass rounding errors
assert _solved.rhs.evalf(7) == (_radial_constant**2 + _axial_constant**2).evalf(7)

# There is no calculate() method. Neutron flux is usually being used internally to pass to other laws.
