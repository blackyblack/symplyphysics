from sympy import Eq, solve, dsolve, Derivative
from symplyphysics import (
    Symbol,
    Function,
    units,
    print_expression,
    Quantity,
    validate_input,
    validate_output,
)

# Description
## The total charge of any isolated system is conserved.

# Law: dq/dt = 0
## q - total charge of the system
## d/dt - derivative with respect to time

# Conditions
## - The system is isolated, i.e. no particles can leave it

total_charge = Function("total_charge", units.charge)
time = Symbol("time", units.time)

law = Eq(Derivative(total_charge(time), time), 0)


def print_law() -> str:
    return print_expression(law)


@validate_input(total_charge_before_=total_charge)
@validate_output(total_charge)
def calculate_charge_after(total_charge_before_: Quantity) -> Quantity:
    dsolved = dsolve(law, total_charge(time))
    dsolved_sub_before = dsolved.subs(total_charge(time), total_charge_before_)
    C1 = solve(dsolved_sub_before, "C1")[0]
    dsolved_sub_after = dsolved.subs("C1", C1)
    result_charge = dsolved_sub_after.rhs
    return Quantity(result_charge)
