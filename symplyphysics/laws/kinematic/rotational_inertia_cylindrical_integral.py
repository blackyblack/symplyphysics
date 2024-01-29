from sympy import Eq, Integral, pi
from symplyphysics import (
    units,
    Symbol,
    Quantity,
    Function,
    validate_input,
    validate_output,
    angle_type,
    print_expression,
)

# Description
## In case of a rigid body with a continuously distributed mass, its rotational inertia is expressed
## as a volume integral over the entire body, i.e. a triple integral over space coordinates.

# Law: I = Integral(rho(V) * r(V)**2, dV) = Integral(rho(r, theta, z) * r**3, dr*dtheta*dz)
## I - rotational inertia
## r - distance to axis of rotation
## rho - density of volume element
## dV = r * dr * dtheta * dz - volume element in cylindrical coordinates (r, theta, z)

# Notes:
## - The integration is carried out over the entire body as to include every volume element.
## - THe z-axis is the rotational axis of the body

rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
density = Function("density", units.mass / units.length**3)

radius = Symbol("radius", units.length)
radius_start = Symbol("radius_start", units.length)
radius_end = Symbol("radius_end", units.length)

polar_angle = Symbol("polar_angle", angle_type)
polar_angle_start = Symbol("polar_angle_start", angle_type)
polar_angle_end = Symbol("polar_angle_end", angle_type)

height = Symbol("height", units.length)
height_start = Symbol("height_start", units.length)
height_end = Symbol("height_end", units.length)

law = Eq(
    rotational_inertia,
    Integral(
        density(radius, polar_angle, height) * radius**3,
        (radius, radius_start, radius_end),
        (polar_angle, polar_angle_start, polar_angle_end),
        (height, height_start, height_end),
    ),
)


def print_law() -> str:
    return print_expression(law)


# Assuming constant density thoughout the body.
# The body is a cylinder with given radius and height.
# The rotational axis is the axis of the cylinder.
@validate_input(density_=density, radius_=radius, height_=height)
@validate_output(rotational_inertia)
def calculate_rotational_inertia(
    density_: Quantity, radius_: Quantity, height_: Quantity
) -> Quantity:
    result = law.rhs.subs({
        density(radius, polar_angle, height): density_,
        radius_start: 0,
        radius_end: radius_,
        polar_angle_start: 0,
        polar_angle_end: 2 * pi,
        height_start: 0,
        height_end: height_,
    }).doit()
    return Quantity(result)
