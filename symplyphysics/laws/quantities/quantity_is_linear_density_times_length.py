"""
Extensive quantity is linear density times length
=================================================

An extensive quantity of interest can be obtained by multiplying the corresponding *linear
density* by length.
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
Extensive quantity. See :symbols:`any_dimension`.
"""

linear_density = clone_as_symbol(
    symbols.any_dimension,
    display_symbol="lambda_X",
    display_latex="\\lambda_X",
)
"""
Intensive linear density, which has the dimension of :attr:`~extensive_quantity`
divided by :attr:`~length`. See :symbols:`any_dimension`.
"""

length = symbols.length
"""
:symbols:`length`.
"""

law = Eq(extensive_quantity, linear_density * length)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(length_=length)
def calculate_extensive_quantity(
    linear_density_: Quantity,
    length_: Quantity,
) -> Quantity:
    result_expr = law.rhs.subs({
        linear_density: linear_density_,
        length: length_,
    })
    result = Quantity(result_expr)

    assert_equivalent_dimension(
        result,
        "return",
        "calculate_extensive_quantity",
        linear_density_.dimension * units.length,
    )

    return result
