from sympy import (I, Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output,)

# Description
## While the serial resistance of ideal capacitor is zero, its impedance depends on its reactance.
## Law: Zc = -j * Xc, where
## Zc is capacitor impedance,
## Xc is capacitive reactance.

capacitor_impedance = Symbol("capacitor_impedance", units.impedance)
capacitive_reactance = Symbol("capacitive_reactance", units.impedance)

law = Eq(capacitor_impedance, -I * capacitive_reactance)


def print_law() -> str:
    return print_expression(law)


@validate_input(reactance_=capacitive_reactance)
@validate_output(capacitor_impedance)
def calculate_impedance(reactance_: Quantity) -> Quantity:
    result_impedance_expr = solve(law, capacitor_impedance, dict=True)[0][capacitor_impedance]
    result_expr = result_impedance_expr.subs({
        capacitive_reactance: reactance_,
    })
    return Quantity(result_expr)
