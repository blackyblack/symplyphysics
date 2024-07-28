"""
Current is charge derivative
============================

The instantaneous electric current, or simply the *electric current*, is a physical quantity
defined as the time derivative of the flowing charge.
"""

from sympy import (Eq, Derivative)
from symplyphysics import (units, Quantity, Function, Symbol, validate_input,
    validate_output)

current = Function("current", units.current)
"""
Electric current as a function of time.

Symbol:
    :code:`I(t)`
"""

charge = Function("charge", units.charge)
"""
Electric charge as a function of time.

Symbol:
    :code:`q(t)`
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

definition = Eq(current(time), Derivative(charge(time), time))
r"""
:code:`I(t) = Derivative(q(t), t)`

Latex:
    .. math::
        I(t) = \frac{d q}{d t}
"""


@validate_input(charge_start_=charge, charge_end_=charge, time_=time)
@validate_output(current)
def calculate_current(charge_start_: Quantity, charge_end_: Quantity, time_: Quantity) -> Quantity:
    charge_function_ = time * (charge_end_ - charge_start_) / time_
    applied_definition = definition.subs(charge(time), charge_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
