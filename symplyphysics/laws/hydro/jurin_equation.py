"""
Capillary rise from surface tension and contact angle
=====================================================

The Jurin's law determines the height to which the liquid rises in capillaries. It
states that the maximum height of a liquid in a capillary tube is inversely proportional
to the tube's diameter.

**Notation:**

#. :quantity_notation:`acceleration_due_to_gravity`.

**Conditions:**

#. The surface of the meniscus is spherical.
#. Height :math:`h` of the raised (lowered) liquid is much larger than the radius
   :math:`r` of the capillary.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Jurin%27s_law>`__.

..
    TODO: rename file to use descriptive name
"""

from sympy import Eq, solve, cos
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities
)

height = symbols.height
"""
:symbols:`height` of the liquid column.
"""

surface_tension = symbols.surface_tension
"""
:symbols:`surface_tension` of the liquid.
"""

angle = symbols.angle
"""
Contact :symbols:`angle` between of the liquid and the tube wall.
"""

density = symbols.density
"""
:symbols:`density` of the liquid.
"""

radius = symbols.radius
"""
:symbols:`radius` of the capillary.
"""

law = Eq(
    height, 2 * surface_tension * cos(angle) /
    (density * radius * quantities.acceleration_due_to_gravity))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(surface_tension_coefficient_=surface_tension,
    angle_=angle,
    density_of_liquid_=density,
    radius_=radius)
@validate_output(height)
def calculate_height(surface_tension_coefficient_: Quantity, angle_: Quantity,
    density_of_liquid_: Quantity, radius_: Quantity) -> Quantity:
    result_expr = solve(law, height, dict=True)[0][height]
    result_height = result_expr.subs({
        surface_tension: surface_tension_coefficient_,
        angle: angle_,
        density: density_of_liquid_,
        radius: radius_,
    })
    result_height_quantity = Quantity(result_height)
    if result_height_quantity.scale_factor < radius_.scale_factor:
        raise ValueError(
            f"The height must be greater than the radius. Currently {result_height_quantity.scale_factor} < {radius_.scale_factor}"
        )
    return result_height_quantity
