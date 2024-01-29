from sympy import Eq, Integral
from symplyphysics import (
    units,
    Symbol,
    Quantity,
    Function,
    validate_input,
    validate_output,
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

rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
density = Function("density", units.mass / units.length**3)
distance_to_axis = Function("distance_to_axis", units.length)

x = Symbol("x", units.length)
x_start = Symbol("x_start", units.length)
x_end = Symbol("x_end", units.length)

y = Symbol("y", units.length)
y_start = Symbol("y_start", units.length)
y_end = Symbol("y_end", units.length)

z = Symbol("z", units.length)
z_start = Symbol("z_start", units.length)
z_end = Symbol("z_end", units.length)

law = Eq(
    rotational_inertia,
    Integral(
    density(x, y, z) * distance_to_axis(x, y, z),
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
        distance_to_axis(x, y, z): x**2 + y**2,
        x_start: -1.0 * x_,
        x_end: x_,
        y_start: -1.0 * y_,
        y_end: y_,
        z_start: -1.0 * z_,
        z_end: z_,
    }).doit()
    return Quantity(result)
