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

heat = Symbol('heat', units.energy)
delta_inner_energy = Symbol('delta_inner_energy', units.energy)
work = Symbol('work', units.energy)

law = Eq(heat, delta_inner_energy + work)


def print_law() -> str:
    return print_expression(law)


@validate_input(delta_inner_energy_=delta_inner_energy, work_=work)
@validate_output(heat)
def calculate_heat(delta_inner_energy_: Quantity, work_: Quantity) -> Quantity:
    """Calculating the law for heat, where the system does the work."""
    solved = solve(law, heat, dict=True)[0][heat]
    result_expr = solved.subs({
        delta_inner_energy: delta_inner_energy_,
        work: work_
    })
    return Quantity(result_expr)
