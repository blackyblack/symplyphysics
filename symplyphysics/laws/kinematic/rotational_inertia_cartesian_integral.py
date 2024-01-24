from sympy import Integral
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
## as an integral over the entire body.

# Law: I = Integral(rho(x, y, z) * (x**2 + y**2), dx*dy*dz)
## I - rotational inertia
## r - perpendicular distance from rotation axis to volume element
## rho - density of the volume element
## x**2 + y**2 = r**2 - square of perpendicular distance to rotation axis
## dx*dy*dz = dV - volume element

# Conditions:
## 1. The z axis is the rotational axis of the body
## 2. The integration is carried out over the entire body as to include every volume element.

density = Function("density", units.mass / units.length**3)
x = Symbol("x", units.length)
y = Symbol("y", units.length)
z = Symbol("z", units.length)


def rotational_inertia_law(x_limits, y_limits, z_limits):
    return Integral(
        density(x, y, z) * (x**2 + y**2),
        (x, *x_limits),
        (y, *y_limits),
        (z, *z_limits),
    )


# Assuming constant density throughout the body
# The body is a cube with sizes 2*x_, 2*y_, 2*z_
@validate_input(density_=density, x_=x, y_=y, z_=z)
@validate_output(units.mass * units.length**2)
def calculate_rotational_inertia(density_: Quantity, x_: Quantity, y_: Quantity, z_: Quantity) -> Quantity:
    result = rotational_inertia_law((-x_, x_), (-y_, y_), (-z_, z_)).subs(
        density(x, y, z), density_
    ).doit()
    return Quantity(result)
