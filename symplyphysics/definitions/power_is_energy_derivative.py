"""
Power is energy derivative
==========================

Power is the amount of energy transferred or converted per unit time. Equally, it is the
rate at which work is done.
"""

from sympy import Eq, Derivative
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
)

power = clone_as_function(symbols.power, display_symbol="P(t)")
"""
:symbols:`power` as a function of time.
"""

energy = clone_as_function(symbols.energy, display_symbol="E(t)")
"""
:symbols:`energy` as a function of time.
"""

time = symbols.time
"""
:symbols:`time`.
"""

definition = Eq(power(time), Derivative(energy(time), time))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(energy_start_=energy, energy_end_=energy, time_=time)
@validate_output(power)
def calculate_power(energy_start_: Quantity, energy_end_: Quantity, time_: Quantity) -> Quantity:
    energy_function_ = time * (energy_end_ - energy_start_) / time_
    applied_definition = definition.subs(energy(time), energy_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
