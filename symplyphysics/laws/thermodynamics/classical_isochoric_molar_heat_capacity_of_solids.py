r"""
Classical isochoric molar heat capacity of solids
=================================================

The *Dulong-Petit law* states that the classical expression for the molar specific heat capacity
of certain chemical elements is constant for temperatures far from the absolute zero.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.

**Conditions:**

#. The temperature of the system is big enough to disregard quantum effects.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Dulong%E2%80%93Petit_law>`__.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_output,
    quantities,
    symbols,
    clone_as_symbol,
)

isochoric_molar_heat_capacity = clone_as_symbol(
    symbols.molar_heat_capacity,
    display_symbol="c_pm",
    display_latex="c_{p, m}",
)
"""
:symbols:`molar_heat_capacity` at constant pressure.
"""

law = Eq(isochoric_molar_heat_capacity, 3 * quantities.molar_gas_constant)
"""
:laws:symbol::

:laws:latex::
"""


@validate_output(isochoric_molar_heat_capacity)
def calculate_isochoric_molar_heat_capacity() -> Quantity:
    return Quantity(law.rhs)
