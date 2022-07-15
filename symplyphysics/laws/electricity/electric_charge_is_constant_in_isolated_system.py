from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Q_after = Q_before
## Where Q_after - summary electric charge of isolated electric system after any period of time, and Q_before - initial summary charge.
## In other words, with no external charge flowing thougth the bouds of system the summary electric charge of all system components remains constant
## during any current flows (charge interchanges) between components. Charge is neither being created in isolated system nor dissapearing.

charge_before, charge_after = symbols('charge_before charge_after')
law = Eq(charge_after, charge_before)

def print():
    return pretty(law, use_unicode=False)

@validate_input(charge_before_ = units.charge)
@validate_output(units.charge)
def calculate_charge_after(charge_before_: Quantity) -> Quantity:
    solved = solve(law, charge_after, dict=True)[0][charge_after]    
    result_expr = solved.subs(charge_before, charge_before_)
    return expr_to_quantity(result_expr, 'charge_after')
