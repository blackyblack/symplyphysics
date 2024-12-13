from sympy import Eq, pi, cos
from sympy.vector import CoordSys3D
from symplyphysics import (Function, Quantity, Symbol, print_expression, units)
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

# TODO: find link

neutron_flux_power_constant = Symbol("C1", 1 / units.area / units.time, constant=True)
x_distance_from_center = Symbol("x_distance_from_center", units.length)
y_distance_from_center = Symbol("y_distance_from_center", units.length)
z_distance_from_center = Symbol("z_distance_from_center", units.length)
parallelepiped_width = Symbol("parallelepiped_width", units.length)
parallelepiped_length = Symbol("parallelepiped_length", units.length)
parallelepiped_height = Symbol("parallelepiped_height", units.length)
neutron_flux = Function("neutron_flux", 1 / units.area / units.time)

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
    neutron_flux(x_distance_from_center, y_distance_from_center, z_distance_from_center),
    neutron_flux_power_constant * cos(width_constant * x_distance_from_center) *
    cos(length_constant * y_distance_from_center) * cos(height_constant * z_distance_from_center))

# Check the solution by passing the known neutron flux to the geometric_buckling_from_neutron_flux.
# Neutron flux is a function of x, y, z in the cartesian coordinates.

# Boundary conditions:
# - vacuum boundary condition: Ф(a + d) = Ф(ae) = Ф(b + d) = Ф(be) = Ф(c + d) = Ф(ce) = 0
# - finite flux condition: 0 <= Ф(x, y, z) < ∞
# - interface condition: the neutron flux and the normal component of the neutron current must be continuous
# - source condition: all neutrons flowing through the bounding area of the source must come from the neutron source

# define flux function in cylindrical coordinates as a function of cylinder radius and height
cartesian_coordinates = CoordSys3D("cartesian_coordinates")
# Make linter happy
x = getattr(cartesian_coordinates, "x")
y = getattr(cartesian_coordinates, "y")
z = getattr(cartesian_coordinates, "z")
unit_length = Quantity(1, dimension=units.length)
neutron_flux_function_cartesian = law.subs({
    x_distance_from_center: x * unit_length,
    y_distance_from_center: y * unit_length,
    z_distance_from_center: z * unit_length
})

solved = geometric_buckling_from_neutron_flux.apply_neutron_flux_function(
    neutron_flux_function_cartesian.rhs)

# check with the derived law: Bg^2 = width_constant**2 + length_constant**2 + height_constant**2
assert solved.rhs == (width_constant**2 + length_constant**2 + height_constant**2)


def print_law() -> str:
    return print_expression(law)


# There is no calculate() method. Neutron flux is usually being used internally to pass to other laws.
