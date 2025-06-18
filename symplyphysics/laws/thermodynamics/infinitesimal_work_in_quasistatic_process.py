"""
Infinitesimal work in quasistatic process
=========================================

For a process in a closed system, occuring slowly enough for it to be quasi-static (i.e. the system
remains in internal physical, but not necessarily chemical, thermodynamic equilibrium), the infinitesimal
increment of work done *by* the system is related to the pressure inside the system and the infinitesimal
increment of volume of the system.

**Notation:**

#. :math:`\\delta` (:code:`delta`) denotes that the increment is an inexact, path-dependent differential.
#. :math:`d` denotes that the increment is an exact, path-independent differential.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Work_(thermodynamics)#Pressure%E2%80%93volume_work>`__.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.core.operations.symbolic import ExactDifferential, InexactDifferential

work_done_by_system = InexactDifferential(symbols.work)
"""
Infinitesimal increment of :symbols:`work` done *by* the system.
"""

pressure = symbols.pressure
"""
:symbols:`pressure` inside the system.
"""

volume_change = ExactDifferential(symbols.volume)
"""
Infinitesimal increment of :symbols:`volume` of the system.
"""

law = Eq(work_done_by_system, pressure * volume_change)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    pressure_inside_system_=pressure,
    infinitesimal_volume_change_=volume_change,
)
@validate_output(work_done_by_system)
def calculate_work_done_by_system(
    pressure_inside_system_: Quantity,
    infinitesimal_volume_change_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        pressure: pressure_inside_system_,
        volume_change: infinitesimal_volume_change_,
    })
    return Quantity(result)
