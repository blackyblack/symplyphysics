from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
                           validate_output)

# Description
# _______________________________________________________
# LAW: Total energy of an isolated system is constant.
# Energy can be transformed from one form to another, but
# can neither be created nor destroyed.
# Formula: dU = Q + A
# Where:
# dU is the change in internal energy of the system,
# A is the work done by the system,
# Q is the heat supplied to the system.
# If the work is done by the system, and not by external
# forces, then the formula is used:
# Q=dU+A
# where A is gas work.

heat = Symbol("heat", units.energy)
delta_inner_energy = Symbol("delta_inner_energy", units.energy)
work = Symbol("work", units.energy)

law = Eq(delta_inner_energy, heat + work)


def print_law() -> str:
    return print_expression(law)


@validate_input(heat_=heat, work_=work)
@validate_output(delta_inner_energy)
def calculate_inner_energy(heat_: Quantity, work_: Quantity) -> Quantity:
    """Calculation of the law for a system where the work
     is performed by external forces. """
    solved = solve(law, delta_inner_energy, dict=True)[0][delta_inner_energy]
    result_expr = solved.subs({
        heat: heat_,
        work: work_
    })
    return Quantity(result_expr)
