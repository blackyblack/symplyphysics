"""
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
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

internal_energy_change = Symbol("internal_energy_change", units.energy)
"""
Infinitesimal change in internal energy of the system.

Symbol:
    :code:`dU`
"""

heat_supplied_to_system = Symbol("heat_supplied_to_system", units.energy)
r"""
Infinitesimal amount of heat supplied to the system during its interaction with the environment.

Symbol:
    :code:`delta Q`

Latex:
    :math:`\delta Q`
"""

work_done_by_system = Symbol("work_done_by_system", units.energy)
r"""
Infinitesimal work done by the system on its environment.

Symbol:
    :code:`delta W`

Latex:
    :math:`\delta W`
"""

law = Eq(internal_energy_change, heat_supplied_to_system - work_done_by_system)
r"""
:code:`dU = delta Q - delta W`

Latex:
    .. math::
        d U = \delta Q - \delta W
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
