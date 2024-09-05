"""
Quantity is volumetric density times volume
===========================================

An extensive quantity of interest can be obtained by multiplying the corresponding *volumetric
density* by volume.
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
Extensive quantity.

Symbol:
    :code:`X`
"""

volumetric_density = SymSymbol("volumetric_density")
r"""
Intensive volumetric density.

Symbol:
    :code:`rho_X`

Latex:
    :math:`\rho_X`
"""

volume = Symbol("volume", units.volume)
"""
Volume.

Symbol:
    :code:`V`
"""

law = Eq(extensive_quantity, volumetric_density * volume)
r"""
:code:`X = rho_X * V`

Latex:
    .. math::
        X = \rho_X V
"""


@validate_input(volume_=volume)
def calculate_extensive_quantity(
    volumetric_density_: Quantity,
    volume_: Quantity,
) -> Quantity:
    result_expr = law.rhs.subs({
        volumetric_density: volumetric_density_,
        volume: volume_,
    })
    result = Quantity(result_expr)

    assert_equivalent_dimension(
        result,
        "return",
        "calculate_extensive_quantity",
        volumetric_density_.dimension * units.volume,
    )

    return result
