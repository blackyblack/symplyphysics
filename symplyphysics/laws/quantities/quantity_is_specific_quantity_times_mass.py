from sympy import Eq, Symbol as SymSymbol
from symplyphysics import (
    units,
    Quantity,
    print_expression,
    validate_input,
    symbols,
)
from symplyphysics.core.dimensions import assert_equivalent_dimension

# Description
## An extensive quantity of interest can be obtained by multiplying the corresponding mass-specific
## intensive quantity by the mass.

# Law: X = x * m
## X - extensive quantity
## x - (intensive) mass-specific quantity
## m - mass

extensive_quantity = SymSymbol("extensive_quantity")
specific_quantity = SymSymbol("specific_quantity")
mass = symbols.basic.mass

law = Eq(extensive_quantity, specific_quantity * mass)


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_=mass)
def calculate_extensive_quantity(
    specific_quantity_: Quantity,
    mass_: Quantity,
) -> Quantity:
    result_expr = law.rhs.subs({
        specific_quantity: specific_quantity_,
        mass: mass_,
    })
    result = Quantity(result_expr)

    assert_equivalent_dimension(
        result,
        "return",
        "calculate_extensive_quantity",
        specific_quantity_.dimension * units.mass,
    )

    return result
