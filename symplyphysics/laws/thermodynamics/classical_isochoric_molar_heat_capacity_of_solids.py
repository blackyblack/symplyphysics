r"""
Classical isochoric molar heat capacity of solids
=================================================

The *Dulong-Petit law* states that the classical expression for the molar specific heat capacity
of certain chemical elements is constant for temperatures far from the absolute zero.

**Notation:**

#. :math:`R` is the molar gas constant.

**Conditions:**

#. The temperature of the system is big enough to disregard quantum effects.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_output,
)

isochoric_molar_heat_capacity = Symbol(
    "isochoric_molar_heat_capacity", units.energy / (units.temperature * units.amount_of_substance))
r"""
Heat capacity at constant volume per unit amount of substance.

Symbol:
    :code:`C_V`

Latex:
    :math:`C_V`
"""

law = Eq(isochoric_molar_heat_capacity, 3 * units.molar_gas_constant)
r"""
:code:`C_V = 3 * R`

Latex:
    .. math::
        C_V = 3 R
"""


@validate_output(isochoric_molar_heat_capacity)
def calculate_isochoric_molar_heat_capacity() -> Quantity:
    return Quantity(law.rhs)
