"""
Current is charge derivative
============================

*Electric current* is the rate at which electric charge changes with time in classical electromagnetism.

**Conditions:**

#. Charge is differentiable with respect to time.

**Links:**

#. `Wikipedia – Electric current <https://en.wikipedia.org/wiki/Electric_current#Metals>`__
"""

from sympy import Eq, Derivative
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

current = clone_as_function(symbols.current, [time])
"""
:symbols:`current` as a function of time.
"""

charge = clone_as_function(symbols.charge, [time])
"""
:symbols:`charge` as a function of time.
"""

definition = Eq(current(time), Derivative(charge(time), time))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(charge_start_=charge, charge_end_=charge, time_=time)
@validate_output(current)
def calculate_current(charge_start_: Quantity, charge_end_: Quantity, time_: Quantity) -> Quantity:
    # charge changes linearly from over the time interval `time_`
    charge_function_ = time * (charge_end_ - charge_start_) / time_
    applied_definition = definition.subs(charge(time), charge_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
