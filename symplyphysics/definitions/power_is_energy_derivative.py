r"""
Power is energy derivative
==========================

Power is the amount of energy transferred or converted per unit time. Equally, it is the
rate at which work is done.

**Notation:**

#. :math:`\frac{d}{d t}` denotes a derivative w.r.t. time.
"""

from sympy import (Eq, Derivative)
from symplyphysics import (units, Quantity, Function, Symbol, validate_input,
    validate_output)

power = Function("power", units.power)
"""
Power as a function of time.

Symbol:
    :code:`P`
"""

energy = Function("energy", units.energy)
"""
Energy as a function of time.

Symbol:
    :code:`E`
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

definition = Eq(power(time), Derivative(energy(time), time))
r"""
:code:`P = dE/dt`

Latex:
    .. math::
        P = \frac{d E}{d t}
"""


@validate_input(energy_start_=energy, energy_end_=energy, time_=time)
@validate_output(power)
def calculate_power(energy_start_: Quantity, energy_end_: Quantity, time_: Quantity) -> Quantity:
    energy_function_ = time * (energy_end_ - energy_start_) / time_
    applied_definition = definition.subs(energy(time), energy_function_)
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
