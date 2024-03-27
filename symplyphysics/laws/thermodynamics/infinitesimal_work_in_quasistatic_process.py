from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## For a process in a closed system, occuring slowly enough for it to be quasi-static (i.e. the system
## remains in internal physical, but not necessarily chemical, thermodynamic equilibrium), the infinitesimal
## increment of work done _by_ the system is related to the pressure inside the system and the infinitesimal
## increment of volume of the system.

# Law: delta W = P * dV
## delta W - an infinitesimal increment of work done _by_ the system
## P - pressure inside the system
## dV - exact differential of an infinitesimal increment of volume of the system

# Note:
## - `delta` means that the increment is an inexact differential.

infinitesimal_work_done = Symbol("infinitesimal_work_done", units.energy)
pressure_inside_system = Symbol("pressure_inside_system", units.pressure)
infinitesimal_volume_change = Symbol("infinitesimal_volume_change", units.volume)

law = Eq(infinitesimal_work_done, pressure_inside_system * infinitesimal_volume_change)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    pressure_inside_system_=pressure_inside_system,
    infinitesimal_volume_change_=infinitesimal_volume_change,
)
@validate_output(infinitesimal_work_done)
def calculate_infinitesimal_work_done(
    pressure_inside_system_: Quantity,
    infinitesimal_volume_change_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        pressure_inside_system: pressure_inside_system_,
        infinitesimal_volume_change: infinitesimal_volume_change_,
    })
    return Quantity(result)
