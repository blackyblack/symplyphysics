"""
Quantity is areal density times volume
======================================

An extensive quantity of interest can be obtained by multiplying the corresponding *areal
density* by area.
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

areal_density = clone_as_symbol(
    symbols.any_dimension,
    display_symbol="sigma_X",
    display_latex="\\sigma_X",
)
r"""
Intensive area-specific density, which has the dimension of :attr:`~extensive_quantity`
divided by :attr:`~area`. See :symbols:`any_dimension`.
"""

area = symbols.area
"""
:symbols:`area`.
"""

law = Eq(extensive_quantity, areal_density * area)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(area_=area)
def calculate_extensive_quantity(
    areal_density_: Quantity,
    area_: Quantity,
) -> Quantity:
    result_expr = law.rhs.subs({areal_density: areal_density_, area: area_})
    result = Quantity(result_expr)

    assert_equivalent_dimension(
        result,
        "return",
        "calculate_extensive_quantity",
        areal_density_.dimension * units.area,
    )

    return result
