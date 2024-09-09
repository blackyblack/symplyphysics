"""
Quantity is molar quantity times amount of substance
====================================================

An extensive quantity of interest can be obtained by multiplying the corresponding molar-specific
intensive quantity by the amount of substance.
"""

from sympy import Eq, Symbol as SymSymbol
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
)
from symplyphysics.core.dimensions import assert_equivalent_dimension

extensive_quantity = SymSymbol("extensive_quantity")
"""
Extensive property.

Symbol:
    :code:`X`
"""

molar_quantity = SymSymbol("molar_quantity")
r"""
Intensive molar quantity.

Symbol:
    :code:`X_m`

Latex:
    :math:`X_m`
"""

amount_of_substance = Symbol("amount_of_substance", units.amount_of_substance)
"""
Amount of substance.

Symbol:
    :code:`n`
"""

law = Eq(extensive_quantity, molar_quantity * amount_of_substance)
r"""
:code:`X = X_m * n`

Latex:
    .. math::
        X = X_m n
"""


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

    assert_equivalent_dimension(
        result,
        "return",
        "calculate_extensive_quantity",
        molar_quantity_.dimension * units.amount_of_substance,
    )

    return result
