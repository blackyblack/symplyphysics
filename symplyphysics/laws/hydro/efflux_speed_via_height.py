"""
Efflux speed via height
=======================

The speed of a fluid flowing out from a small orifice can be expressed as a function of
the height of the fluid column. It is also known as the **Torricelli's law**.

**Notation:**

#. :quantity_notation:`acceleration_due_to_gravity`.

**Conditions:**

#. The orifice is very small relative to the horizontal cross-section of the container.
#. The fluid is :ref:`ideal <ideal_fluid_def>`.
#. The fluid is subjected to the gravity force of Earth.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Torricelli%27s_law#>`__.
"""

from sympy import Eq, solve, sqrt
from symplyphysics import Quantity, validate_input, validate_output, symbols, quantities

efflux_speed = symbols.flow_speed
"""
:symbols:`flow_speed` of the fluid flowing out of the pipe.
"""

height = symbols.height
"""
:symbols:`height` of the fluid column above the orifice.
"""

law = Eq(efflux_speed, sqrt(2 * quantities.acceleration_due_to_gravity * height))
"""
:laws:symbol::

:laws:latex::
"""

# TODO Derive from Bernoulli's equation and constancy of volume flux


@validate_input(height_=height)
@validate_output(efflux_speed)
def calculate_velocity(height_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, efflux_speed, dict=True)[0][efflux_speed]
    result_expr = result_velocity_expr.subs({height: height_})
    return Quantity(result_expr)
