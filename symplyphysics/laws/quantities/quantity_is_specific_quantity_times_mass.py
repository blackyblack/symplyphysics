"""
Quantity is specific quantity times mass
========================================

An extensive quantity of interest can be obtained by multiplying the corresponding mass-specific
intensive quantity by the mass.
"""

from sympy import Eq, Symbol as SymSymbol
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    symbols,
)
from symplyphysics.core.dimensions import assert_equivalent_dimension

extensive_quantity = SymSymbol("extensive_quantity")
"""
Extensive quantity.

Symbol:
    :code:`X`
"""

specific_quantity = SymSymbol("specific_quantity")
r"""
Intensive mass-specific quantity.

Symbol:
    :code:`x`
"""

mass = symbols.mass
"""
:symbols:`mass`.

Symbol:
    :code:`m`
"""

law = Eq(extensive_quantity, specific_quantity * mass)
r"""
:code:`X = x * m`

Latex:
    .. math::
        X = x m
"""


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
