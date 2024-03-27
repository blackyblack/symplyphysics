from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
    # symbols,
)

# Description
## The first law of thermodynamics is a formulation of the law of conservation of energy in the
## context of thermodynamic processes. It relates the change in internal energe of the system,
## work done on or by the system and the amount of heat supplied to or widthrawn from the system
## during a thermodynamic process.

# Law: dU = delta Q - delta W
## dU - infinitesimal change in internal energy of system
## delta Q - inexact differential of an infinitesimal amount of heat supplied to
##           or withdrawn from the system during the interaction with its surroundings
## delta W - inexact differential of infinitesimal work done by the system

# Note
## - The sign before work `W` is negative for work done _by_ the system and positive for 
##   work done  _on_ the system by the environment.
## - The sign before the amount of heat `Q` is positive for heat flowing from the environment 
##   _into_ the system and negative for heat flowing _out of_ the system into the environment.
## - This formula can be extended to finite processes, with the only difference being that the
##   changes are not infinitesimal but finite (`dU` would refer to a finite change in the system's
##   internal energy, delta Q would be the total amount of heat supplied and delta W would be the total
##   work done by the system).

# Conditions
## - This formula applies to any infinitesimal thermodynamic process, be it (quasi)static or not.

internal_energy_change = Symbol("internal_energy_change", units.energy)
heat_supplied_to_system = Symbol("heat_supplied_to_system", units.energy)
work_done_by_system = Symbol("work_done_by_system", units.energy)

law = Eq(internal_energy_change, heat_supplied_to_system - work_done_by_system)


def print_law() -> str:
    return print_expression(law)


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
