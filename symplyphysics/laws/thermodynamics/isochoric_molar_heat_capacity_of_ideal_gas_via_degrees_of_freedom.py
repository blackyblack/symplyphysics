r"""
Isochoric molar heat capacity of ideal gas via degrees of freedom
=================================================================

The internal energy of an ideal gas consisting of molecules whose energy is all kinetic
depends only on the degrees of freedom of the molecule and the temperature of the gas.
From that one can derive the expression of the isochoric heat capacity of ideal gases.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.

**Notes:**

#. For applications, see :doc:`internal energy of ideal gas <laws.thermodynamics.internal_energy_of_ideal_gas_via_temperature>`.
#. :math:`f = 3` for monatomic molecules.
#. :math:`f = 5` for diatomic molecules.
#. :math:`f = 6` for non-linear polyatomic molecules.

**Conditions:**

#. Gas is ideal.
#. Works in the classical theory of heat capacity of gases. For a more accurate representation refer to
   the quantum theory, which accounts for the "freezing" of the degrees of freedom and other phenomena.

**Links:**

#. `Physics LibreTexts, "Equipartition Theorem" <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_II_-_Thermodynamics_Electricity_and_Magnetism_(OpenStax)/02%3A_The_Kinetic_Theory_of_Gases/2.04%3A_Heat_Capacity_and_Equipartition_of_Energy>`__.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    quantities,
    symbols,
    clone_as_symbol,
)

isochoric_molar_heat_capacity = clone_as_symbol(
    symbols.molar_heat_capacity,
    display_symbol="c_Vm",
    display_latex="c_{V, \\text{m}}",
)
"""
:symbols:`molar_heat_capacity` at constant :symbols:`volume`.
"""

degrees_of_freedom = symbols.degrees_of_freedom
"""
Number of :symbols:`degrees_of_freedom` of gas molecules.
"""

law = Eq(isochoric_molar_heat_capacity, (degrees_of_freedom / 2) * quantities.molar_gas_constant)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(degrees_of_freedom_=degrees_of_freedom)
@validate_output(isochoric_molar_heat_capacity)
def calculate_isochoric_molar_heat_capacity(degrees_of_freedom_: int) -> Quantity:
    result = law.rhs.subs(degrees_of_freedom, degrees_of_freedom_)
    return Quantity(result)
