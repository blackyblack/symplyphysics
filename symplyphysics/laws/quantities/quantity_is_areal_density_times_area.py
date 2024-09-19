"""
Quantity is areal density times volume
======================================

An extensive quantity of interest can be obtained by multiplying the corresponding *areal
density* by area.
"""

from sympy import Eq, Symbol as SymSymbol
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
)

extensive_quantity = SymSymbol("extensive_quantity")
"""
Extensive quantity.

Symbol:
    :code:`X`
"""

areal_density = SymSymbol("areal_density")
r"""
Intensive area-specific density.

Symbol:
    :code:`sigma_X`

Latex:
    :math:`\sigma_X`
"""

area = SymbolNew("A", units.area)
"""
Area.
"""

law = Eq(extensive_quantity, areal_density * area)
r"""
:code:`X = sigma_X * A`

Latex:
    .. math::
        X = \sigma_X A
"""

@validate_input(area_=area)
def calculate_extensive_quantity(
    areal_density_: Quantity,
    area_: Quantity,
) -> Quantity:
    result_expr = law.rhs.subs({
        areal_density: areal_density_,
        area: area_
    })
    result = Quantity(result_expr)

    return result
