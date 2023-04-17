from sympy.vector import CoordSys3D
from symplyphysics import (
    symbols, Function, Eq, pretty, sin, pi
)
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_neutron_flux

# Description
## Neutron flux formula for the uniform spherical reactor. The spherical reactor is situated in spherical
## geometry at the origin of coordinates.

## Law: Ф(r) = C1 * (sin(Pi * r / R) / r)
## Where:
## C1 - neutron flux power constant.
## r - distance from the center of sphere.
## Pi - Pi constant.
## R - sphere radius.
## Ф(r) - neutron flux density.

neutron_flux_power_constant = symbols("C1", constant=True)
distance_from_center = symbols("distance_from_center")
sphere_radius = symbols("sphere_radius")
neutron_flux_function = symbols("neutron_flux_function", cls = Function)

# This constant is being used for geometric buckling calculation
# See: [geometric buckling for uniform sphere](geometric_buckling_for_uniform_sphere.py)
radial_constant = pi / sphere_radius

law = Eq(neutron_flux_function(distance_from_center),
    neutron_flux_power_constant * sin(radial_constant * distance_from_center) / distance_from_center)

# Check the solution by passing the known neutron flux to the geometric_buckling_from_neutron_flux.
# Neutron flux is a function of radius in the spherical coordinates.

# Boundary conditions:
# - vacuum boundary condition: Ф(R + d) = Ф(Re) = 0
# - finite flux condition: 0 <= Ф(r) < ∞
# - interface condition: the neutron flux and the normal component of the neutron current must be continuous
# - source condition: all neutrons flowing through the bounding area of the source must come from the neutron source
# - albedo boundary condition: Ф(Ralbedo) = 0

# define flux function in spherical coordinates as a function of sphere radius
spherical_coordinates = CoordSys3D("spherical_coordinates", transformation="spherical")
neutron_flux_function_spherical = law.subs(distance_from_center, spherical_coordinates.r)

solved = geometric_buckling_from_neutron_flux.apply_neutron_flux_function(neutron_flux_function_spherical.rhs)

# check with the derived law: Bg^2 = radial_constant**2
assert solved.rhs == radial_constant**2

def print():
    return pretty(law, use_unicode=False)

# There is no calculate() method. Neutron flux is usually being used internally to pass to other laws.
