from sympy import Eq, Symbol as SymSymbol
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
)

# Description
## An extensive quantity of interest can be obtained by multiplying the corresponding molar-specific
## intensive quantity by the amount of substance.

# Law: X = X_m * n
## X - extensive quantity
## X_m - (intensive) molar quantity
## n - amount of substance

extensive_quantity = SymSymbol("extensive_quantity")
molar_quantity = SymSymbol("molar_quantity")
amount_of_substance = Symbol("amount_of_substance", units.amount_of_substance)

law = Eq(extensive_quantity, molar_quantity * amount_of_substance)


def print_law() -> str:
    return print_expression(law)


@validate_input(amount_of_substance_=amount_of_substance)
def calculate_extensive_quantity(
    molar_quantity_: Quantity,
    amount_of_substance_: Quantity,
) -> Quantity:
    result_expr = law.rhs.subs({
        molar_quantity: molar_quantity_,
        amount_of_substance: amount_of_substance_,
    })
    result = Quantity(result_expr)

    return result
