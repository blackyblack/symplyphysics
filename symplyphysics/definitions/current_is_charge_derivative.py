"""
Current is charge derivative
============================

The instantaneous electric current, or simply the *electric current*, is a physical quantity
defined as the time derivative of the flowing charge.
"""

from sympy import (Eq, Derivative)
from symplyphysics import (units, Quantity, FunctionNew, SymbolNew, validate_input, validate_output)

current = FunctionNew("I(t)", units.current, display_latex="I")
"""
Electric current as a function of time.
"""

charge = FunctionNew("q(t)", units.charge, display_latex="q")
"""
Electric charge as a function of time.
"""

time = SymbolNew("t", units.time)
"""
Time.
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
