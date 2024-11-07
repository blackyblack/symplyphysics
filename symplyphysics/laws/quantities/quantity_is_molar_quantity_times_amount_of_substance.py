"""
Quantity is molar quantity times amount of substance
====================================================

An extensive quantity of interest can be obtained by multiplying the corresponding molar-specific
intensive quantity by the amount of substance.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.dimensions import assert_equivalent_dimension

extensive_quantity = symbols.any_dimension
"""
Extensive property. See :symbols:`any_dimension`.
"""

molar_quantity = clone_as_symbol(
    symbols.any_dimension,
    display_symbol="x_m",
    display_latex="x_m",
)
"""
Intensive molar quantity, which has the dimension of :attr:`~extensive_quantity`
divided by :attr:`~amount_of_substance`. See :symbols:`any_dimension`.
"""

amount_of_substance = symbols.amount_of_substance
"""
:symbols:`amount_of_substance`.
"""

law = Eq(extensive_quantity, molar_quantity * amount_of_substance)
"""
:laws:symbol::

:laws:latex::
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
