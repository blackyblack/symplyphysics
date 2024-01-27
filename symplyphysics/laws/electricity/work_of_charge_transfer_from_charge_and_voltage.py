from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## When the test charge moves in an electric field, we can talk about the work being done
## at the moment by electric forces. This operation depends on the amount of charge and voltage.

## Law is: A = q * U, where
## A - work of charge transfer,
## q - electric charge,
## U - voltage.

work_of_charge_transfer = Symbol("work_of_charge_transfer", units.energy)

charge = Symbol("charge", units.charge)
voltage = Symbol("voltage", units.voltage)

law = Eq(work_of_charge_transfer, charge * voltage)


def print_law() -> str:
    return print_expression(law)


@validate_input(charge_=charge, voltage_=voltage)
@validate_output(work_of_charge_transfer)
def calculate_work(charge_: Quantity, voltage_: Quantity) -> Quantity:
    result_expr = solve(law, work_of_charge_transfer, dict=True)[0][work_of_charge_transfer]
    result_expr = result_expr.subs({
        charge: charge_,
        voltage: voltage_,
    })
    return Quantity(result_expr)
