"""
Fractional change is change over initial value
==============================================


"""

from sympy import Eq, Symbol as SymSymbol
from symplyphysics import (
    dimensionless,
    Quantity,
    Symbol,
    validate_output,
)
from symplyphysics.core.dimensions import assert_equivalent_dimension

fractional_change = Symbol("fractional_change", dimensionless)
r"""
Fractional change of quantity :math:`X`

Symbol:
    :code:`e_X`

Latex:
    :math:`e_X`
"""

change = SymSymbol("change")
r"""
Change in the value of the quantity.

Symbol:
    :code:`dX`

Latex:
    :math:`\Delta X`
"""

initial_value = SymSymbol("initial_value")
"""
Initial value of the quantity.

Symbol:
    :code:`X`
"""

law = Eq(fractional_change, change / initial_value)
r"""
:code:`e_X = dX / X`

Latex:
    .. math::
        e_X = \frac{\Delta X}{X}
"""


@validate_output(fractional_change)
def calculate_fractional_change(
    change_: Quantity,
    initial_value_: Quantity,
) -> Quantity:
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
    
    return Quantity(result)
