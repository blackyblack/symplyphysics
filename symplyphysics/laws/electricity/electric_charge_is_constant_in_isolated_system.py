from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## Q_after = Q_before
## Where Q_after - summary electric charge of isolated electric system after any period of time, and Q_before - initial summary charge.
## In other words, with no external charge flowing thougth the bouds of system the summary electric charge of all system components remains constant
## during any current flows (charge interchanges) between components. Charge is neither being created in isolated system nor dissapearing.

charge_before = Symbol("charge_before", units.charge)
charge_after = Symbol("charge_after", units.charge)

law = Eq(charge_after, charge_before)


def print_law() -> str:
    return print_expression(law)


@validate_input(charge_before_=charge_before)
@validate_output(charge_after)
def calculate_charge_after(charge_before_: Quantity) -> Quantity:
    solved = solve(law, charge_after, dict=True)[0][charge_after]
    result_expr = solved.subs(charge_before, charge_before_)
    return Quantity(result_expr)
