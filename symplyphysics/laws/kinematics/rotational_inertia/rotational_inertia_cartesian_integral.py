"""
Rotational inertia in terms of Cartesian integral
=================================================

In case of a rigid body with a continuously distributed mass, its rotational inertia is expressed
as a volume integral over the entire body, i.e. a triple integral over :math:`x, y, z` in Cartesian
coordinates.

**Notes:**

#. The integration is carried out over the entire body as to include every volume element.

**Links:**

#. `Wikipedia, derivable from fourth equation <https://en.wikipedia.org/wiki/Moment_of_inertia#Point_mass>`__.
"""

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

rotational_inertia = symbols.rotational_inertia
"""
:symbols:`rotational_inertia`.
"""

x = symbols.position
"""
:symbols:`position` on the :math:`x` axis.
"""

x_start = clone_as_symbol(x, subscript="0")
"""
Initial position on the :math:`x` axis.
"""

x_end = clone_as_symbol(x, subscript="1")
"""
Final position on the :math:`x` axis.
"""

y = clone_as_symbol(symbols.position, display_symbol="y", display_latex="y")
"""
:symbols:`position` on the :math:`y` axis.
"""

y_start = clone_as_symbol(y, subscript="0")
"""
Initial position on the :math:`y` axis.
"""

y_end = clone_as_symbol(y, subscript="1")
"""
Final position on the :math:`y` axis.
"""

z = clone_as_symbol(symbols.position, display_symbol="z", display_latex="z")
"""
:symbols:`position` on the :math:`z` axis.
"""

z_start = clone_as_symbol(z, subscript="0")
"""
Initial position on the :math:`z` axis.
"""

z_end = clone_as_symbol(z, subscript="1")
"""
Final position on the :math:`z` axis.
"""

density = clone_as_function(symbols.density, [x, y, z])
"""
Mass-specific :symbols:`density` as a function of :attr:`~x`, :attr:`~y`, :attr:`~z`.
"""

distance_to_axis = clone_as_function(symbols.distance_to_axis, [x, y, z])
"""
:symbols:`distance_to_axis` as a function of :attr:`~x`, :attr:`~y`, :attr:`~z`.
"""

law = Eq(
    rotational_inertia,
    Integral(
    density(x, y, z) * distance_to_axis(x, y, z)**2,
    (x, x_start, x_end),
    (y, y_start, y_end),
    (z, z_start, z_end),
    ),
)
"""
:laws:symbol::

:laws:latex::
"""


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


# UNIQUE_LAW_ID: 465
