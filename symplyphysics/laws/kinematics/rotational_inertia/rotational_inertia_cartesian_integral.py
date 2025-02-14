from sympy import Eq, Integral, sqrt
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
    clone_as_function,
)

# Description
## In case of a rigid body with a continuously distributed mass, its rotational inertia is expressed
## as a volume integral over the entire body, i.e. a triple integral over x, y, and z in Cartesian
## coordinates.

# Law: I = Integral(rho(V) * r(V)**2, dV) = Integral(rho(x,y,z) * r(x,y,z)**2, dx*dy*dz)
## I - rotational inertia
## r - perpendicular distance from rotation axis to volume element
## rho - density of the volume element
## dV = dx*dy*dz - volume element, where (x, y, z) are Cartesian coordinates

# Note:
## - The integration is carried out over the entire body as to include every volume element.

# Links:
## Wikipedia, derivable from fourth equation <https://en.wikipedia.org/wiki/Moment_of_inertia#Point_mass>
# TODO: update documentation

rotational_inertia = symbols.rotational_inertia

x = symbols.position
x_start = clone_as_symbol(x, subscript="0")
x_end = clone_as_symbol(x, subscript="1")

y = clone_as_symbol(symbols.position, display_symbol="y", display_latex="y")
y_start = clone_as_symbol(y, subscript="0")
y_end = clone_as_symbol(y, subscript="1")

z = clone_as_symbol(symbols.position, display_symbol="z", display_latex="z")
z_start = clone_as_symbol(z, subscript="0")
z_end = clone_as_symbol(z, subscript="1")

density = clone_as_function(symbols.density, [x, y, z])
distance_to_axis = clone_as_function(symbols.distance_to_axis, [x, y, z])

law = Eq(
    rotational_inertia,
    Integral(
    density(x, y, z) * distance_to_axis(x, y, z)**2,
    (x, x_start, x_end),
    (y, y_start, y_end),
    (z, z_start, z_end),
    ),
)


# Assuming constant density throughout the body.
# The body is a rectangular parallelepiped with sizes 2*x_, 2*y_, 2*z_.
# The rotational axis passes through its center and is parallel to the z-axis.
@validate_input(density_=density, x_=x, y_=y, z_=z)
@validate_output(units.mass * units.length**2)
def calculate_rotational_inertia(density_: Quantity, x_: Quantity, y_: Quantity,
    z_: Quantity) -> Quantity:
    result = law.rhs.subs({
        density(x, y, z): density_,
        distance_to_axis(x, y, z): sqrt(x**2 + y**2),
        x_start: -1.0 * x_,
        x_end: x_,
        y_start: -1.0 * y_,
        y_end: y_,
        z_start: -1.0 * z_,
        z_end: z_,
    }).doit()
    return Quantity(result)
