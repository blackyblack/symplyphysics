"""
Inner pressure is sum of pressures
==================================

Inner pressure of an ideal fluid is the sum of static, dynamic, and hydrostatic
pressures at a chosen point in space.

**Conditions:**

#. The fluid is :ref:`ideal <ideal_fluid_def>`.

**Links:**

#. `Wikipedia, derivable from here <https://en.wikipedia.org/wiki/Bernoulli%27s_principle#Incompressible_flow_equation>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

inner_pressure = Symbol("inner_pressure", units.pressure)
r"""
Inner pressure of the fluid.

Symbol:
    :code:`p_inner`

Latex:
    :math:`p_\text{inner}`
"""

static_pressure = Symbol("static_pressure", units.pressure)
r"""
Static pressure of the fluid.

Symbol:
    :code:`p_static`

Latex:
    :math:`p_\text{static}`
"""

dynamic_pressure = Symbol("dynamic_pressure", units.pressure)
"""
Dynamic pressure of the fluid.

Symbol:
    :code:`q`
"""

hydrostatic_pressure = Symbol("hydrostatic_pressure", units.pressure)
"""
Hydrostatic pressure of the fluid.

Symbol:
    :code:`p`
"""

law = Eq(inner_pressure, static_pressure + dynamic_pressure + hydrostatic_pressure)
r"""
:code:`p_inner = p_static + q + p`

Latex:
    .. math::
        p_\text{inner} = p_\text{static} + q + p
"""


@validate_input(
    static_pressure_=static_pressure,
    dynamic_pressure_=dynamic_pressure,
    hydrostatic_pressure_=hydrostatic_pressure,
)
@validate_output(inner_pressure)
def calculate_inner_pressure(
    static_pressure_: Quantity,
    dynamic_pressure_: Quantity,
    hydrostatic_pressure_: Quantity,
) -> Quantity:
    result_expr = solve(law, inner_pressure)[0]
    result_inner_pressure = result_expr.subs({
        static_pressure: static_pressure_,
        dynamic_pressure: dynamic_pressure_,
        hydrostatic_pressure: hydrostatic_pressure_,
    })
    return Quantity(result_inner_pressure)
