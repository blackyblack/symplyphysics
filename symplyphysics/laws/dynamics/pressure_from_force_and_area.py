"""
Pressure from force and area
============================

Pressure is the amount of force applied perpendicular to the surface of an object
per unit area over which the force is distributed.
"""

from sympy import (Eq, solve)
from symplyphysics import (symbols, units, Quantity, Symbol, validate_input, validate_output)

pressure = Symbol("pressure", units.pressure)
"""
Pressure.

Symbol:
    :code:`p`
"""

force = symbols.force
"""
:attr:`~symplyphysics.symbols.force` applied.
"""

area = Symbol("area", units.area)
r"""
Area over which the force is distributed.

Symbol: 
    :code:`A`
"""

law = Eq(pressure, force / area)
r"""
:code:`p = F / A`

Latex:
    .. math::
        p = \frac{F}{A}
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
