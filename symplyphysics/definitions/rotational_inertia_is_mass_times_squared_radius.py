"""
Rotational inertia is mass times squared radius
===============================================

Rotational inertia, or moment of inertia, is a physical quantity that describes a body's ability to
be inertial during rotation. It is a rotational analog of mass in linear motion.

**Conditions:**

#. The object is a material point, rigid and uniform.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, symbols)

rotational_inertia = Symbol("rotational_inertia", units.mass * units.area)
"""
Rotational inertia of the material point.

Symbol:
    :code:`I`
"""

mass = symbols.basic.mass
"""
:attr:`~symplyphysics.symbols.basic.mass` of the material point.
"""

radial_distance = Symbol("radial_distance", units.length)
"""
Distance to the axis of rotation.

Symbol:
    :code:`r`
"""

definition = Eq(rotational_inertia, mass * radial_distance**2)
"""
:code:`I = m * r^2`

Latex:
    .. math::
        I = m r^2
"""


@validate_input(mass_=mass, radius_=radial_distance)
@validate_output(rotational_inertia)
def calculate_moment_of_inertia(mass_: Quantity, radius_: Quantity) -> Quantity:
    result_inertia_expr = solve(definition, rotational_inertia, dict=True)[0][rotational_inertia]
    result_expr = result_inertia_expr.subs({mass: mass_, radial_distance: radius_})
    return Quantity(result_expr)
