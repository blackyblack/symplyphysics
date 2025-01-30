r"""
Internal energy change via heat and work
========================================

The *first law of thermodynamics* is a formulation of the law of conservation of energy in the
context of thermodynamic processes. It relates the change in internal energe of the system,
work done on or by the system and the amount of heat supplied to or withdrawn from the system
during a thermodynamic process.

**Notation:**

#. :math:`\delta` (:code:`delta`) denotes an inexact, path-dependent differential.
#. :math:`d` denotes an exact, path-independent differential.

**Notes:**

#. This formula applies to any infinitesimal thermodynamic process, be in static or not.
#. This formula can be extended to finite processes by changing infinitesimal changes to finite ones.
#. The work done *by* the system onto the environment is positive, and work done *onto* the system
   by the environment is negative.
#. The heat flowing from the environment *into* the system is positive, and the heat flowing into
   the environment *out of* the system is negative.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/First_law_of_thermodynamics#Definition>`__.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

internal_energy_change = clone_as_symbol(symbols.internal_energy, display_symbol="dU", display_latex="dU")
"""
Infinitesimal change in :symbols:`internal_energy` of the system.
"""

heat_supplied_to_system = clone_as_symbol(symbols.heat, display_symbol="delta(Q)", display_latex="\\delta Q")
r"""
Infinitesimal amount of :symbols:`heat` supplied to the system during its interaction
with the environment.
"""

work_done_by_system = clone_as_symbol(symbols.work, display_symbol="delta(W)", display_latex="\\delta W")
"""
Infinitesimal :symbols:`work` done by the system on its environment.
"""

law = Eq(internal_energy_change, heat_supplied_to_system - work_done_by_system)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    heat_supplied_to_system_=heat_supplied_to_system,
    work_done_by_system_=work_done_by_system,
)
@validate_output(internal_energy_change)
def calculate_internal_energy_change(
    heat_supplied_to_system_: Quantity,
    work_done_by_system_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        heat_supplied_to_system: heat_supplied_to_system_,
        work_done_by_system: work_done_by_system_,
    })
    return Quantity(result)
