r"""
Infinitesimal work in quasistatic process
=========================================

For a process in a closed system, occuring slowly enough for it to be quasi-static (i.e. the system
remains in internal physical, but not necessarily chemical, thermodynamic equilibrium), the infinitesimal
increment of work done *by* the system is related to the pressure inside the system and the infinitesimal
increment of volume of the system.

**Notation:**

#. :math:`\delta` (:code:`delta`) denotes that the increment is an inexact, path-dependent differential.
#. :math:`d` denotes that the increment is an exact, path-independent differential.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Work_(thermodynamics)#Pressure%E2%80%93volume_work>`__.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

infinitesimal_work_done = Symbol("infinitesimal_work_done", units.energy)
r"""
Infinitesimal increment of work done *by* the system.

Symbol:
    :code:`delta(W)`

Latex:
    :math:`\delta W`
"""

pressure = Symbol("pressure", units.pressure)
"""
Pressure inside the system.

Symbol:
    :code:`p`
"""

infinitesimal_volume_change = Symbol("infinitesimal_volume_change", units.volume)
"""
Infinitesimal increment of volume of the system.

Symbol:
    :code:`dV`

Latex:
    :math:`dV`
"""

law = Eq(infinitesimal_work_done, pressure * infinitesimal_volume_change)
r"""
:code:`delta(W) = p * dV`

Latex:
    .. math::
        \delta W = p \, dV
"""


@validate_input(
    pressure_inside_system_=pressure,
    infinitesimal_volume_change_=infinitesimal_volume_change,
)
@validate_output(infinitesimal_work_done)
def calculate_infinitesimal_work_done(
    pressure_inside_system_: Quantity,
    infinitesimal_volume_change_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        pressure: pressure_inside_system_,
        infinitesimal_volume_change: infinitesimal_volume_change_,
    })
    return Quantity(result)
