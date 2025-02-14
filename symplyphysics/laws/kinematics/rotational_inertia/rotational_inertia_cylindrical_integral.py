from sympy import Eq, Integral, pi
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
    clone_as_function,
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
## - The z-axis is the rotational axis of the body

# Links:
## Wikipedia, derivable from fourth equation <https://en.wikipedia.org/wiki/Moment_of_inertia#Point_mass>
# TODO: update documentation

rotational_inertia = symbols.rotational_inertia

radius = clone_as_symbol(symbols.radius)
radius_start = clone_as_symbol(radius, subscript="0")
radius_end = clone_as_symbol(radius, subscript="1")

polar_angle = clone_as_symbol(symbols.angle)
polar_angle_start = clone_as_symbol(polar_angle, subscript="0")
polar_angle_end = clone_as_symbol(polar_angle, subscript="1")

height = clone_as_symbol(symbols.height)
height_start = clone_as_symbol(height, subscript="0")
height_end = clone_as_symbol(height, subscript="1")

density = clone_as_function(symbols.density, [radius, polar_angle, height])

law = Eq(
    rotational_inertia,
    Integral(
    density(radius, polar_angle, height) * radius**3,
    (radius, radius_start, radius_end),
    (polar_angle, polar_angle_start, polar_angle_end),
    (height, height_start, height_end),
    ),
)


# Assuming constant density throughout the body.
# The body is a cylinder with given radius and height.
# The rotational axis is the axis of the cylinder.
@validate_input(density_=density, radius_=radius, height_=height)
@validate_output(rotational_inertia)
def calculate_cylinder_rotational_inertia(density_: Quantity, radius_: Quantity,
    height_: Quantity) -> Quantity:
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
