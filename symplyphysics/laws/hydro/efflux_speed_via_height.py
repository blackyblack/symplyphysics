"""
Efflux speed via height
=======================

The speed of a fluid flowing out from a small orifice can be expressed as a function
of the height of the fluid column. It is also known as the *Torricelli's law*.

**Conditions:**

#. The orifice is very small relative to the horizontal cross-section of the container.
#. The fluid is :ref:`ideal <ideal_fluid_def>`.
#. The fluid is subjected to the gravity force of Earth.
"""

from sympy import (Eq, solve, sqrt)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

efflux_speed = Symbol("efflux_speed", units.velocity)
"""
Speed of the fluid flowing out of the pipe.

Symbol:
    :code:`v`
"""

height = Symbol("height", units.length)
"""
Height of the fluid column above the orifice.

Symbol:
    :code:`h`
"""

law = Eq(efflux_speed, sqrt(2 * units.acceleration_due_to_gravity * height))
r"""
:code:`v = sqrt(2 * g * h)`

Latex:
    .. math::
        v = \sqrt{2 g h}
"""

# TODO Derive from Bernoulli's equation and constancy of volume flux


@validate_input(height_=height)
@validate_output(efflux_speed)
def calculate_velocity(height_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, efflux_speed, dict=True)[0][efflux_speed]
    result_expr = result_velocity_expr.subs({height: height_})
    return Quantity(result_expr)
