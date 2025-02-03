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
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

inner_pressure = clone_as_symbol(symbols.pressure, display_symbol="p_inner", display_latex="p_\\text{inner}")
"""
Inner :symbols:`pressure` of the fluid.
"""

static_pressure = clone_as_symbol(symbols.pressure, display_symbol="p_static", display_latex="p_\\text{static}")
"""
Static :symbols:`pressure` of the fluid.
"""

dynamic_pressure = symbols.dynamic_pressure
"""
:symbols:`dynamic_pressure` of the fluid.
"""

hydrostatic_pressure = symbols.hydrostatic_pressure
"""
:symbols:`hydrostatic_pressure` of the fluid.
"""

law = Eq(inner_pressure, static_pressure + dynamic_pressure + hydrostatic_pressure)
"""
:laws:symbol::

:laws:latex::
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
