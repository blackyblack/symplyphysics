"""
Current is charge derivative
============================

The instantaneous electric current, or simply the *electric current*, is a physical quantity
defined as the time derivative of the flowing charge.
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
Time.
"""

current = clone_as_function(symbols.current, [time])
"""
Electric :symbols:`current` as a function of time.
"""

charge = clone_as_function(symbols.charge, [time])
"""
Electric :symbols:`charge` as a function of time.
"""

definition = Eq(current(time), Derivative(charge(time), time))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(charge_start_=charge, charge_end_=charge, time_=time)
@validate_output(current)
def calculate_current(charge_start_: Quantity, charge_end_: Quantity, time_: Quantity) -> Quantity:
    charge_function_ = time * (charge_end_ - charge_start_) / time_
    applied_definition = definition.subs(charge(time), charge_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
