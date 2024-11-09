"""
Quantity is volumetric density times volume
===========================================

An extensive quantity of interest can be obtained by multiplying the corresponding *volumetric
density* by volume.
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

extensive_quantity = symbols.any_quantity
"""
Extensive quantity. See :symbols:`any_dimension`.
"""

volumetric_density = clone_as_symbol(
    symbols.any_quantity,
    display_symbol="rho_X",
    display_latex="\\rho_X",
)
"""
Intensive volumetric density, which has the dimension of :attr:`~extensive_quantity`
divided by :attr:`~volume`. See :symbols:`any_dimension`.
"""

volume = symbols.volume
"""
:symbols:`volume`.
"""

law = Eq(extensive_quantity, volumetric_density * volume)
"""
:laws:symbol::

:laws:latex::
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
