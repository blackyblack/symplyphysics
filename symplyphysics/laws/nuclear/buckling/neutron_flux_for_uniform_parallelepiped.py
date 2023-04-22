from sympy import Eq, symbols, pi, Function as SymFunction, cos
from sympy.vector import CoordSys3D
from symplyphysics import print_expression
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_neutron_flux
from symplyphysics.laws.nuclear.buckling import neutron_flux_for_uniform_slab

# Description
## Neutron flux formula for the uniform rectangular parallelepiped reactor of physical dimensions a * b * c.
## The parallelepiped reactor is situated in cartesian geometry at the origin of coordinates.

## Law: Ф(x, y, z) = C1 * cos(Pi * x / a) * cos(Pi * y / b) * cos(Pi * z / c)
## Where:
## C1 - neutron flux power constant.
## Pi - Pi constant.
## x, y, z - distances (width, length, height) from the center of parallelepiped.
## a - width of the parallelepiped.
## b - length of the parallelepiped.
## c - height of the parallelepiped.
## Ф(x, y, z) - neutron flux density.

neutron_flux_power_constant = symbols("C1", constant=True)
x_distance_from_center = symbols("x_distance_from_center")
y_distance_from_center = symbols("y_distance_from_center")
z_distance_from_center = symbols("z_distance_from_center")
parallelepiped_width = symbols("parallelepiped_width")
parallelepiped_length = symbols("parallelepiped_length")
parallelepiped_height = symbols("parallelepiped_height")
neutron_flux_function = symbols("neutron_flux_function", cls=SymFunction)

# These constants are being used for geometric buckling calculation
# See: [geometric buckling for uniform parallelepiped](geometric_buckling_for_uniform_parallelepiped.py)
width_constant = pi / parallelepiped_width
length_constant = pi / parallelepiped_length
height_constant = pi / parallelepiped_height

# derived the same way as uniform slab axial_constant
assert width_constant == neutron_flux_for_uniform_slab.axial_constant.subs(
    neutron_flux_for_uniform_slab.slab_width, parallelepiped_width)
assert length_constant == neutron_flux_for_uniform_slab.axial_constant.subs(
    neutron_flux_for_uniform_slab.slab_width, parallelepiped_length)
assert height_constant == neutron_flux_for_uniform_slab.axial_constant.subs(
    neutron_flux_for_uniform_slab.slab_width, parallelepiped_height)

law = Eq(
    neutron_flux_function(x_distance_from_center, y_distance_from_center,
                          z_distance_from_center),
    neutron_flux_power_constant * cos(width_constant * x_distance_from_center) *
    cos(length_constant * y_distance_from_center) *
    cos(height_constant * z_distance_from_center))

# Check the solution by passing the known neutron flux to the geometric_buckling_from_neutron_flux.
# Neutron flux is a function of x, y, z in the cartesian coordinates.

# Boundary conditions:
# - vacuum boundary condition: Ф(a + d) = Ф(ae) = Ф(b + d) = Ф(be) = Ф(c + d) = Ф(ce) = 0
# - finite flux condition: 0 <= Ф(x, y, z) < ∞
# - interface condition: the neutron flux and the normal component of the neutron current must be continuous
# - source condition: all neutrons flowing through the bounding area of the source must come from the neutron source

# define flux function in cylindrical coordinates as a function of cylinder radius and height
cartesian_coordinates = CoordSys3D("cartesian_coordinates")
neutron_flux_function_cartesian = law.subs({
    x_distance_from_center: cartesian_coordinates.x,
    y_distance_from_center: cartesian_coordinates.y,
    z_distance_from_center: cartesian_coordinates.z
})

solved = geometric_buckling_from_neutron_flux.apply_neutron_flux_function(
    neutron_flux_function_cartesian.rhs)

# check with the derived law: Bg^2 = width_constant**2 + length_constant**2 + height_constant**2
assert solved.rhs == (width_constant**2 + length_constant**2 +
                      height_constant**2)


def print() -> str:
    return print_expression(law)


# There is no calculate() method. Neutron flux is usually being used internally to pass to other laws.
