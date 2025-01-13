from sympy import (Eq, solve, pi)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output)
from symplyphysics.laws.nuclear.buckling import neutron_flux_for_uniform_cylinder as cylinder_flux

# Description
## Geometric buckling for the uniform cylindrical reactor of physical radius R and height H.
## The cylindrical reactor is situated in cylindrical geometry at the origin of coordinates.

## Law: Bg^2cylinder = (2.405 / R)^2 + (Pi / H)^2
## Where:
## Pi - Pi constant.
## R - radius of the cylinder.
## H - height of the cylinder.
## Bg^2cylinder - squared geometric buckling for cylinder.
##   See [geometric buckling](./geometric_buckling_from_neutron_flux.py) implementation.

# Links:
## Wikipedia, second row in first table <https://en.wikipedia.org/wiki/Geometric_and_material_buckling#Geometric_Buckling>

cylinder_radius = Symbol("cylinder_radius", units.length)
cylinder_height = Symbol("cylinder_height", units.length)
geometric_buckling_squared = Symbol("geometric_buckling_squared", 1 / units.area)

law = Eq(geometric_buckling_squared, (2.405 / cylinder_radius)**2 + (pi / cylinder_height)**2)

# This law is derived from geometric buckling definition (see geometric_buckling_from_neutron_flux.py),
# neutron flux laplacian in cylindrical coordinates and boundary conditions.

# Unfortunately sympy does not support solving with complex boundary conditions so we simply check with known
# solution for the neutron flux:
# See [neutron flux for uniform cylinder](./neutron_flux_for_uniform_cylinder.py)
geometric_buckling_cylinder_squared = cylinder_flux.radial_constant**2 + cylinder_flux.axial_constant**2
geometric_buckling_cylinder_solved = geometric_buckling_cylinder_squared.subs({
    cylinder_flux.cylinder_radius: cylinder_radius,
    cylinder_flux.cylinder_height: cylinder_height
})
assert geometric_buckling_cylinder_solved.evalf(7) == law.rhs.evalf(7)


@validate_input(cylinder_radius_=cylinder_radius, cylinder_height_=cylinder_height)
@validate_output(geometric_buckling_squared)
def calculate_geometric_buckling_squared(cylinder_radius_: Quantity,
    cylinder_height_: Quantity) -> Quantity:
    solved = solve(law, geometric_buckling_squared, dict=True)[0][geometric_buckling_squared]
    result_expr = solved.subs({
        cylinder_radius: cylinder_radius_,
        cylinder_height: cylinder_height_
    })
    return Quantity(result_expr)
