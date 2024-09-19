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
#. Works in the classical theory of heat capacity of gases. For a more accurate represention refer to
   the quantum theory, which accounts for the "freezing" of the degrees of freedom and other phenomena.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    quantities,
)

isochoric_molar_heat_capacity = Symbol(
    "isochoric_molar_heat_capacity",
    units.energy / (units.temperature * units.amount_of_substance),
)
r"""
Heat capacity of ideal gas at constant volume per unit amount of substance.

Symbol:
    :code:`C_V`

Latex:
    :math:`C_V`
"""

degrees_of_freedom = Symbol("degrees_of_freedom", integer=True)
"""
Number of degrees of freedom of gas molecules.

Symbol:
    :code:`f`
"""

law = Eq(isochoric_molar_heat_capacity, (degrees_of_freedom / 2) * quantities.molar_gas_constant)
r"""
:code:`C_V = (f / 2) * R`

Latex:
    .. math::
        C_V = \frac{f}{2} R
"""


@validate_input(degrees_of_freedom_=degrees_of_freedom)
@validate_output(isochoric_molar_heat_capacity)
def calculate_isochoric_molar_heat_capacity(degrees_of_freedom_: int) -> Quantity:
    result = law.rhs.subs(degrees_of_freedom, degrees_of_freedom_)
    return Quantity(result)
