from sympy import sin
from sympy.vector import CoordSys3D
from symplyphysics import (
    symbols, Function, Eq, pretty
)
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_neutron_flux

# Description
## Neutron flux formula for the uniform spherical reactor. The spherical reactor is situated in spherical
## geometry at the origin of coordinates.

## Law: Ф(r) = C1 * (sin(Bgsphere * r) / r)
## Where:
## C1 - neutron flux power constant.
## r - distance from the center of sphere.
## Bgsphere - geometric buckling for sphere.
##   See [geometric buckling](./geometric_buckling_from_neutron_flux.py) implementation.
## Ф(x) - neutron flux density.

neutron_flux_power_constant = symbols('C1', constant=True)
distance_from_center = symbols('distance_from_center')
geometric_buckling_sphere = symbols('geometric_buckling_sphere')
neutron_flux_function = symbols('neutron_flux_function', cls = Function)

law = Eq(neutron_flux_function(distance_from_center),
    neutron_flux_power_constant * sin(geometric_buckling_sphere * distance_from_center) / distance_from_center)

# Check the solution by passing the known neutron flux to the geometric_buckling_from_neutron_flux.
# Neutron flux is a function of radius in the spherical coordinates.

# define flux function in spherical coordinates as a function of sphere radius
spherical_coordinates = CoordSys3D('spherical_coordinates', transformation='spherical')
neutron_flux_function_spherical = law.subs(distance_from_center, spherical_coordinates.r)

solved = geometric_buckling_from_neutron_flux.calculate_geometric_buckling_squared(neutron_flux_function_spherical.rhs)

# check with the derived law
assert solved.rhs == geometric_buckling_sphere**2

def print():
    return pretty(law, use_unicode=False)

# There is no calculate() method. Neutron flux is usually being used internally to pass to the other laws.
