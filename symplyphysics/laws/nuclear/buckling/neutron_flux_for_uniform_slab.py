from sympy import cos, pi
from sympy.vector import CoordSys3D
from symplyphysics import (
    symbols, Function, Eq, pretty
)
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_neutron_flux

# Description
## Neutron flux formula for the reactor in the shape of a slab of physical width a in the x-direction
## and infinite in the y- and z-directions.

## Law: Ф(x) = C1 * cos(Pi * x / a)
## Where:
## C1 - neutron flux power constant.
## x - distance from the center of slab.
## Pi - Pi constant.
## a - slab width.
## Ф(x) - neutron flux density.

neutron_flux_power_constant = symbols('C1', constant=True)
distance_from_center = symbols('distance_from_center')
slab_width = symbols('slab_width')
neutron_flux_function = symbols('neutron_flux_function', cls = Function)

# This constant is being used for geometric buckling calculation
# See: [geometric buckling for uniform slab](geometric_buckling_for_uniform_slab.py)
axial_constant = pi / slab_width

law = Eq(neutron_flux_function(distance_from_center),
    neutron_flux_power_constant * cos(axial_constant * distance_from_center))

# Check the solution by passing the known neutron flux to the geometric_buckling_from_neutron_flux.
# Neutron flux is a function of x coordinate in the cartesian coordinates.

# Boundary conditions:
# - vacuum boundary condition: Ф(a / 2 + d) = Ф(ae / 2) = 0
# - finite flux condition: 0 <= Ф(x) < ∞
# - interface condition: the neutron flux and the normal component of the neutron current must be continuous
# - source condition: all neutrons flowing through the bounding area of the source must come from the neutron source

# define flux function in cartesian coordinates as a function of x coordinate
cartesian_coordinates = CoordSys3D('cartesian_coordinates')
neutron_flux_function_cartesian = law.subs(distance_from_center, cartesian_coordinates.x)

solved = geometric_buckling_from_neutron_flux.apply_neutron_flux_function(neutron_flux_function_cartesian.rhs)

# check with the derived law: Bg^2 = axial_constant**2
assert solved.rhs == axial_constant**2

def print():
    return pretty(law, use_unicode=False)

# There is no calculate() method. Neutron flux is usually being used internally to pass to the other laws.
