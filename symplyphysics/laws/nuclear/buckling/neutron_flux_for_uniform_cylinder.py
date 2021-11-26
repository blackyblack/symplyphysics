from sympy import cos, pi
from sympy.functions.special.bessel import besselj
from sympy.vector import CoordSys3D
from symplyphysics import (
    symbols, Function, Eq, pretty
)
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_neutron_flux

# Description
## Neutron flux formula for the uniform cylindrical reactor of physical radius R and height H.
## The cylindrical reactor is situated in cylindrical geometry at the origin of coordinates.

## Law: Ф(r, z) = C1 * J0(2.405 * r / R) * cos(Pi * z / H)
## Where:
## C1 - neutron flux power constant.
## Pi - Pi constant.
## J0 - Bessel function of the first kind, of the 0 order.
## r - radial distance from the center of cylinder.
## z - vertical distance from the center of cylinder.
## R - radius of the cylinder.
## H - height of the cylinder.
## Ф(r, z) - neutron flux density.

neutron_flux_power_constant = symbols('C1', constant=True)
radial_distance_from_center = symbols('radial_distance_from_center')
vertical_distance_from_center = symbols('vertical_distance_from_center')
cylinder_radius, cylinder_height = symbols('cylinder_radius cylinder_height')
neutron_flux_function = symbols('neutron_flux_function', cls = Function)

radial_constant = 2.405 / cylinder_radius
axial_constant = pi / cylinder_height

law = Eq(neutron_flux_function(radial_distance_from_center, vertical_distance_from_center),
    neutron_flux_power_constant * besselj(0, radial_constant * radial_distance_from_center) *
    cos(axial_constant * vertical_distance_from_center))

# Check the solution by passing the known neutron flux to the geometric_buckling_from_neutron_flux.
# Neutron flux is a function of radius and height in the cylindrical coordinates.

# define flux function in cylindrical coordinates as a function of cylinder radius and height
cylindrical_coordinates = CoordSys3D('cylindrical_coordinates', transformation='cylindrical')
neutron_flux_function_cylindrical = law.subs({
    radial_distance_from_center: cylindrical_coordinates.r,
    vertical_distance_from_center: cylindrical_coordinates.z})

solved = geometric_buckling_from_neutron_flux.calculate_geometric_buckling_squared(neutron_flux_function_cylindrical.rhs)

# check with the derived law: Bg^2 = radial_constant**2 + axial_constant**2
# limit decimals to bypass rounding errors
assert solved.rhs.evalf(7) == (radial_constant**2 + axial_constant**2).evalf(7)

def print():
    return pretty(law, use_unicode=False)

# There is no calculate() method. Neutron flux is usually being used internally to pass to the other laws.
