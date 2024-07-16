"""
Extensive quantity is linear density times length
=================================================

An extensive quantity of interest can be obtained by multiplying the corresponding linear
density by length.
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
Extensive quantity

Symbol:
    X
"""

linear_density = SymSymbol("linear_density")
r"""
Intensive linear density

Symbol:
    lambda_X

Latex:
    :math:`\lambda_X`
"""

length = Symbol("length", units.length)
"""
Length

Symbol:
    L
"""

law = Eq(extensive_quantity, linear_density * length)
r"""
X = lambda_X * L

Latex:
    :math:`X = \lambda_X L`
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
