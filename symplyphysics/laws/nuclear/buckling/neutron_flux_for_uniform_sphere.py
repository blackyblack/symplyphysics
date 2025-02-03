"""
Neutron flux for uniform sphere
===============================

Neutron flux for uniform sphere is a function of radial_distance :math:`r` to the center of the
sphere and depends on the radius :math:`r_0` of the sphere.

**Links:**

#. `NuclearPower, see end of page <https://www.nuclear-power.com/nuclear-power/reactor-physics/neutron-diffusion-theory/finite-spherical-reactor/>`__.
"""

from sympy import Eq, pi, sin
from sympy.vector import CoordSys3D
from symplyphysics import Quantity, units, symbols, clone_as_symbol
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_neutron_flux

dimension_factor = clone_as_symbol(symbols.neutron_flux, subscript="0")
"""
Dimension factor that appears as a coefficient in the solution to the :ref:`differential
equation <Diffusion equation from neutron flux>`. See :symbols:`neutron_flux`.
"""

radial_distance = symbols.distance_to_origin
"""
Radial distance, or :symbols:`distance_to_origin` of coordinate system, i.e. the center
of the sphere.
"""

radius = clone_as_symbol(symbols.radius, subscript="0")
"""
:symbols:`radius` of the sphere.
"""

neutron_flux = symbols.neutron_flux
"""
:symbols:`neutron_flux` at a :attr:`radial_distance` from the center of the sphere.
"""

# This constant is being used for geometric buckling calculation
# See: [geometric buckling for uniform sphere](geometric_buckling_for_uniform_sphere.py)
_radial_constant = pi / radius

law = Eq(neutron_flux, dimension_factor * (sin(_radial_constant * radial_distance) / radial_distance))
"""
:laws:symbol::

:laws:latex::
"""

# Check the solution by passing the known neutron flux to the geometric_buckling_from_neutron_flux.
# Neutron flux is a function of radius in the spherical coordinates.

# Boundary conditions:
# - vacuum boundary condition: Ф(R + d) = Ф(Re) = 0
# - finite flux condition: 0 <= Ф(r) < ∞
# - interface condition: the neutron flux and the normal component of the neutron current must be continuous
# - source condition: all neutrons flowing through the bounding area of the source must come from the neutron source
# - albedo boundary condition: Ф(Ralbedo) = 0

# define flux function in spherical coordinates as a function of sphere radius
_spherical_coordinates = CoordSys3D("_spherical_coordinates", transformation="spherical")
# Make linter happy
_r = getattr(_spherical_coordinates, "r")
_unit_length = Quantity(1, dimension=units.length)
_neutron_flux_function_spherical = law.subs(radial_distance, _r * _unit_length)

_solved = geometric_buckling_from_neutron_flux.apply_neutron_flux_function(
    _neutron_flux_function_spherical.rhs)

# check with the derived law: Bg^2 = _radial_constant**2
assert _solved.rhs == _radial_constant**2

# There is no calculate() method. Neutron flux is usually being used internally to pass to other laws.
