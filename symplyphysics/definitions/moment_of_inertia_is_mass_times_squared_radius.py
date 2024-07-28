"""
Moment of inertia is mass times squared radius
==============================================

Moment of inertia, or rotational inertia, is a physical quantity that describes a body's ability to
be inertial during rotation. It is a rotational analog of :attr:`~symplyphysics.symbols.basic.mass` in
linear motion.

**Conditions:**

#. The object is a material point, rigid and uniform.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, symbols)

moment_of_inertia = Symbol("moment_of_inertia", units.mass * units.area)
"""
Moment of inertia of the material point.

Symbol:
    :code:`I`
"""

mass = symbols.basic.mass
"""
:attr:`~symplyphysics.symbols.basic.mass` of the material point.

Symbol:
    :code:`m`
"""

spinning_radius = Symbol("spinning_radius", units.length)
"""
Distance to the axis of rotation.

Symbol:
    :code:`r`
"""

definition = Eq(moment_of_inertia, mass * spinning_radius**2)
"""
:code:`I = m * r^2`

Latex:
    .. math::
        I = m r^2
"""


@validate_input(mass_=symbols.basic.mass, radius_=spinning_radius)
@validate_output(moment_of_inertia)
def calculate_moment_of_inertia(mass_: Quantity, radius_: Quantity) -> Quantity:
    result_inertia_expr = solve(definition, moment_of_inertia, dict=True)[0][moment_of_inertia]
    result_expr = result_inertia_expr.subs({mass: mass_, spinning_radius: radius_})
    return Quantity(result_expr)
