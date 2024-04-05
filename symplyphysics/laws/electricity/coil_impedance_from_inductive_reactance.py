from sympy import (I, Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output,)

# Description
## The impedance of ideal coil depends on its inductive reactance. While having zero resistivity, the real part of
## coil impedance is zero.
## Law: Zl = -j * Xl, where
## Zl is coil impedance,
## Xl is inductive reactance.

coil_impedance = Symbol("coil_impedance", units.impedance)
inductive_reactance = Symbol("inductive_reactance", units.impedance)

law = Eq(coil_impedance, I * inductive_reactance)


def print_law() -> str:
    return print_expression(law)


@validate_input(reactance_=inductive_reactance)
@validate_output(coil_impedance)
def calculate_impedance(reactance_: Quantity) -> Quantity:
    result_impedance_expr = solve(law, coil_impedance, dict=True)[0][coil_impedance]
    result_expr = result_impedance_expr.subs({
        inductive_reactance: reactance_,
    })
    return Quantity(result_expr)
