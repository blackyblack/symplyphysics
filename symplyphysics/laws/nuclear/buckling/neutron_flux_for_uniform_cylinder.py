from sympy.functions.special.bessel import besselj
from sympy.vector import CoordSys3D
from symplyphysics import (
    symbols, Function, Eq, pretty, cos, pi
)
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_neutron_flux
from symplyphysics.laws.nuclear.buckling import neutron_flux_for_uniform_slab

# Description
## Neutron flux formula for the uniform cylindrical reactor of physical radius R and height H.
## The cylindrical reactor is situated in cylindrical geometry at the origin of coordinates.

## Law: Ф(r, z) = C1 * J0(2.405 * r / R) * cos(Pi * z / H)
## Where:
## C1 - neutron flux power constant.
## Pi - Pi constant.
## J0 - Bessel function of the first kind, of the 0 order.
## r - radial distance from the center of cylinder.
## z - axial distance from the center of cylinder.
## R - radius of the cylinder.
## H - height of the cylinder.
## Ф(r, z) - neutron flux density.

neutron_flux_power_constant = symbols("C1", constant=True)
radial_distance_from_center = symbols("radial_distance_from_center")
axial_distance_from_center = symbols("axial_distance_from_center")
cylinder_radius, cylinder_height = symbols("cylinder_radius cylinder_height")
neutron_flux_function = symbols("neutron_flux_function", cls = Function)

# These constants are being used for geometric buckling calculation
# See: [geometric buckling for uniform cylinder](geometric_buckling_for_uniform_cylinder.py)
radial_constant = 2.405 / cylinder_radius
axial_constant = pi / cylinder_height

# derived the same way as uniform slab axial_constant
assert axial_constant == neutron_flux_for_uniform_slab.axial_constant.subs(neutron_flux_for_uniform_slab.slab_width, cylinder_height)

law = Eq(neutron_flux_function(radial_distance_from_center, axial_distance_from_center),
    neutron_flux_power_constant * besselj(0, radial_constant * radial_distance_from_center) *
    cos(axial_constant * axial_distance_from_center))

# Check the solution by passing the known neutron flux to the geometric_buckling_from_neutron_flux.
# Neutron flux is a function of radius and height in the cylindrical coordinates.

# Boundary conditions:
# - vacuum boundary condition: Ф(R + d) = Ф(Re) = Ф(a + d) = Ф(ae) = 0
# - finite flux condition: 0 <= Ф(r, x) < ∞
# - interface condition: the neutron flux and the normal component of the neutron current must be continuous
# - source condition: all neutrons flowing through the bounding area of the source must come from the neutron source
# - albedo boundary condition: Ф(Ralbedo) = 0

# radial_constant is the solution of the Bessel function J0, with a condition the neutron flux cannot
# have negative values (finite flux condition) and with zero flux boundary condition

# define flux function in cylindrical coordinates as a function of cylinder radius and height
cylindrical_coordinates = CoordSys3D("cylindrical_coordinates", transformation="cylindrical")
neutron_flux_function_cylindrical = law.subs({
    radial_distance_from_center: cylindrical_coordinates.r,
    axial_distance_from_center: cylindrical_coordinates.z})

solved = geometric_buckling_from_neutron_flux.apply_neutron_flux_function(neutron_flux_function_cylindrical.rhs)

# check with the derived law: Bg^2 = radial_constant**2 + axial_constant**2
# limit decimals to bypass rounding errors
assert solved.rhs.evalf(7) == (radial_constant**2 + axial_constant**2).evalf(7)

def print():
    return pretty(law, use_unicode=False)

# There is no calculate() method. Neutron flux is usually being used internally to pass to other laws.
