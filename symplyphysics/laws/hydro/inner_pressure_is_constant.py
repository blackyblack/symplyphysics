"""
Inner pressure is constant
==========================

Bernoulli's equation applied to an ideal liquid specifies that the inner
pressure of the fluid is constant at all points along a streamline.

**Conditions:**

#. The fluid must be :ref:`ideal <ideal_fluid_def>`.
"""

from sympy import Eq, dsolve, Derivative
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
)

inner_pressure = Function("inner_pressure", units.pressure)
r"""
Inner presure of the fluid at a chosen point in space.

Symbol:
    :code:`p_inner(t)`

Latex:
    :math:`p_\text{inner}(t)`

..
    TODO add link to definition
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

law = Eq(Derivative(inner_pressure(time), time), 0)
r"""
:code:`Derivative(p_inner(t), t) = 0`

Latex:
    .. math::
        \frac{d p_\text{inner}}{d t} = 0
"""


@validate_input(inner_pressure_before_=inner_pressure)
@validate_output(inner_pressure)
def calculate_inner_pressure(inner_pressure_before_: Quantity) -> Quantity:
    dsolved = dsolve(law, inner_pressure(time))
    result_expr = dsolved.subs("C1", inner_pressure_before_).rhs
    return Quantity(result_expr)