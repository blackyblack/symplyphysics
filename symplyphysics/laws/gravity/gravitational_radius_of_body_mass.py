"""
Gravitational radius of massive body
====================================

The **gravitational radius** is a characteristic radius defined for any physical body with mass.
This is the radius of the sphere on which the event horizon created by this mass would be
located (from the point of view of general theory of relativity) if it were distributed
spherically symmetrically, would be stationary (in particular, it would not rotate, but
radial movements are permissible) and would lie entirely inside this sphere.

**Notation:**

#. :quantity_notation:`gravitational_constant`

#. :quantity_notation:`speed_of_light`

..
    TODO rename file
"""

from sympy import Eq, solve
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
    quantities,
)

gravitational_radius = symbols.radius
"""
Gravitational :symbols:`radius` of the body.
"""

body_mass = symbols.mass
"""
:symbols:`mass` of the body.
"""

law = Eq(gravitational_radius, 2 * quantities.gravitational_constant * body_mass / quantities.speed_of_light**2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(body_mass_=body_mass)
@validate_output(gravitational_radius)
def calculate_radius(body_mass_: Quantity) -> Quantity:
    result_expr = solve(law, gravitational_radius, dict=True)[0][gravitational_radius]
    result_expr = result_expr.subs({
        body_mass: body_mass_,
    })
    return Quantity(result_expr)
