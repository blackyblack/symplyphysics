"""
Charge is constant
==================

The total charge of any isolated system is conserved.

**Conditions:**

#. The system is isolated, i.e. no particles can leave it.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Charge_conservation>`__.
"""

from sympy import Eq, solve, dsolve, Derivative
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    clone_as_function,
    symbols,
)

# Law: dq/dt = 0
## q - total charge of the system
## d/dt - derivative with respect to time

time = symbols.time
"""
:symbols:`time`.
"""

total_charge = clone_as_function(symbols.charge, [time])
"""
:symbols:`charge` as a function of :attr:`~time`.
"""

law = Eq(Derivative(total_charge(time), time), 0)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(total_charge_before_=total_charge)
@validate_output(total_charge)
def calculate_charge_after(total_charge_before_: Quantity) -> Quantity:
    dsolved = dsolve(law, total_charge(time))
    dsolved_sub_before = dsolved.subs(total_charge(time), total_charge_before_)
    C1 = solve(dsolved_sub_before, "C1")[0]
    dsolved_sub_after = dsolved.subs("C1", C1)
    result_charge = dsolved_sub_after.rhs
    return Quantity(result_charge)
