"""
Pressure from force and area
============================

Pressure is the amount of force applied perpendicular to the surface of an object
per unit area over which the force is distributed.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Pressure#Formula>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
)

pressure = symbols.pressure
"""
:symbols:`pressure`.
"""

force = symbols.force
"""
:symbols:`force` applied.
"""

area = symbols.area
r"""
:symbols:`area` over which the force is distributed.
"""

law = Eq(pressure, force / area)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(force_=force, area_=area)
@validate_output(pressure)
def calculate_pressure(force_: Quantity, area_: Quantity) -> Quantity:
    result_expr = solve(law, pressure, dict=True)[0][pressure]
    result_pressure = result_expr.subs({
        force: force_,
        area: area_,
    })

    return Quantity(result_pressure)
