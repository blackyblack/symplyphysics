"""
Fractional change is change over initial value
==============================================

Fractional change of a quantity is defined as ratio of the linear change of the quantity
to the initial value of the quantity.
"""

from sympy import Eq
from symplyphysics import (
    convert_to_float,
    Quantity,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.dimensions import assert_equivalent_dimension

fractional_change = clone_as_symbol(symbols.fractional_change, subscript="X")
r"""
:symbols:`fractional_change` of quantity :math:`X`.
"""

change = clone_as_symbol(
    symbols.any_quantity,
    display_symbol="Delta(X)",
    display_latex="\Delta X"
)
"""
Change in the value of the quantity. See :symbols:`any_dimension`.
"""

initial_value = symbols.any_quantity
"""
Initial value of the quantity. See :symbols:`any_dimension`.
"""

law = Eq(fractional_change, change / initial_value)
"""
:laws:symbol::

:laws:latex::
"""


@validate_output(fractional_change)
def calculate_fractional_change(
    change_: Quantity,
    initial_value_: Quantity,
) -> float:
    assert_equivalent_dimension(
        change_,
        "change_",
        "calculate_fractional_change",
        initial_value_,
    )

    result = law.rhs.subs({
        change: change_,
        initial_value: initial_value_,
    })

    return convert_to_float(result)
