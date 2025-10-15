"""
Rotational inertia in terms of a cylindrical integral
=====================================================

In case of a rigid body with a continuously distributed mass, its rotational inertia is expressed
as a volume integral over the entire body, i.e. a triple integral over space coordinates.

**Notes:**

#. The integration is carried out over the entire body as to include every volume element.

**Conditions:**

#. The :math:`z`-axis is the rotational axis of the body.

**Links:**

#. `Wikipedia, derivable from fourth equation <https://en.wikipedia.org/wiki/Moment_of_inertia#Point_mass>`__.
"""

from sympy import Eq, Integral, pi
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
    clone_as_function,
)

rotational_inertia = symbols.rotational_inertia
"""
:symbols:`rotational_inertia` of the body.
"""

radius = clone_as_symbol(symbols.radius)
"""
:symbols:`radius`, or distance to the rotational axis.
"""

radius_start = clone_as_symbol(radius, subscript="0")
"""
Initial :symbols:`radius`.
"""

radius_end = clone_as_symbol(radius, subscript="1")
"""
Final :symbols:`radius`.
"""

polar_angle = clone_as_symbol(symbols.angle)
"""
Polar :symbols:`angle`.
"""

polar_angle_start = clone_as_symbol(polar_angle, subscript="0")
"""
Initial polar :symbols:`angle`.
"""

polar_angle_end = clone_as_symbol(polar_angle, subscript="1")
"""
Final polar :symbols:`angle`.
"""

height = clone_as_symbol(symbols.height)
"""
:symbols:`height`.
"""

height_start = clone_as_symbol(height, subscript="0")
"""
Initial :symbols:`height`.
"""

height_end = clone_as_symbol(height, subscript="1")
"""
Final :symbols:`height`.
"""

density = clone_as_function(symbols.density, [radius, polar_angle, height])
"""
:symbols:`density` as a function of :attr:`~radius`, :attr:`~polar_angle`, and :attr:`~height`.
"""

law = Eq(
    rotational_inertia,
    Integral(
    density(radius, polar_angle, height) * radius**3,
    (radius, radius_start, radius_end),
    (polar_angle, polar_angle_start, polar_angle_end),
    (height, height_start, height_end),
    ),
)
"""
:laws:symbol::

:laws:latex::
"""


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


# UNIQUE_LAW_ID: 466
