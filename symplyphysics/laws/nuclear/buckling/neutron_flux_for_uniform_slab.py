from sympy import cos
from sympy.vector import CoordSys3D
from symplyphysics import (
    symbols, Function, Eq, pretty
)
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_neutron_flux

# Description
## Neutron flux formula for the reactor in the shape of a slab of physical width a in the x-direction
## and infinite in the y- and z-directions.

## Law: Ф(r) = C1 * cos(Bgslab * x)
## Where:
## C1 - neutron flux power constant.
## x - distance from the center of slab.
## Bgslab - geometric buckling for slab.
##   See [geometric buckling](./geometric_buckling_from_neutron_flux.py) implementation.

neutron_flux_power_constant = symbols('C1', constant=True)
distance_from_center = symbols('distance_from_center')
geometric_buckling_slab = symbols('geometric_buckling_slab')
neutron_flux_function = symbols('neutron_flux_function', cls = Function)

law = Eq(neutron_flux_function(distance_from_center),
    neutron_flux_power_constant * cos(geometric_buckling_slab * distance_from_center))

# Check the solution by passing the known neutron flux to the geometric_buckling_from_neutron_flux.
# Neutron flux is a function of x coordinate in the cartesian coordinates.

# define flux function in cartesian coordinates as a function of x coordinate
cartesian_coordinates = CoordSys3D('cartesian_coordinates')
neutron_flux_function_cartesian = law.subs(distance_from_center, cartesian_coordinates.x)

solved = geometric_buckling_from_neutron_flux.calculate_geometric_buckling_squared(neutron_flux_function_cartesian.rhs)

# check with the derived law
assert solved.rhs == geometric_buckling_slab**2

def print():
    return pretty(law, use_unicode=False)

# There is no calculate() method. Neutron flux is usually being used internally to pass to the other laws.
