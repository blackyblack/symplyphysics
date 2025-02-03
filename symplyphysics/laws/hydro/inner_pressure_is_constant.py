"""
Inner pressure is constant
==========================

Bernoulli's equation applied to an ideal liquid specifies that the inner pressure of the
fluid is constant at all points along a streamline.

**Conditions:**

#. The fluid must be :ref:`ideal <ideal_fluid_def>`.

**Links:**

#. `Wikipedia, derivable from here <https://en.wikipedia.org/wiki/Bernoulli%27s_principle#Simplified_form>`__.
"""

from sympy import Eq, dsolve, Derivative
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
)

time = symbols.time
"""
:symbols:`time`.
"""

inner_pressure = clone_as_function(symbols.pressure, [time], display_symbol="p_inner", display_latex="p_\\text{inner}")
"""
Inner pressure of the fluid at a chosen point in space as a function of :attr:`~time`.
See :ref:`Inner pressure is sum of pressures`.
"""

law = Eq(Derivative(inner_pressure(time), time), 0)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(inner_pressure_before_=inner_pressure)
@validate_output(inner_pressure)
def calculate_inner_pressure(inner_pressure_before_: Quantity) -> Quantity:
    dsolved = dsolve(law, inner_pressure(time))
    result_expr = dsolved.subs("C1", inner_pressure_before_).rhs
    return Quantity(result_expr)
