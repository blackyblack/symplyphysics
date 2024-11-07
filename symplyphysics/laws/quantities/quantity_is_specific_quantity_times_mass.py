"""
Quantity is specific quantity times mass
========================================

An extensive quantity of interest can be obtained by multiplying the corresponding mass-specific
intensive quantity by the mass.
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

specific_quantity = clone_as_symbol(
    symbols.any_dimension,
    display_symbol="x",
    display_latex="x",
)
"""
Intensive mass-specific quantity, which has the dimension of :attr:`~extensive_quantity`
divided by :attr:`~mass`. See :symbols:`any_dimension`.
"""

mass = symbols.mass
"""
:symbols:`mass`.
"""

law = Eq(extensive_quantity, specific_quantity * mass)
"""
:laws:symbol::

:laws:latex::
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
