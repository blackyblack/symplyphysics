r"""
Pressure is Maclaurin series of strain
======================================

A statement more general than the Hooke's law is that pressure in a body can be expressed
as a *power series* (Maclaurin series, to be exact) of (engineering normal) *strain* with the
free coefficient being zero, since pressure disappears with the disappearance of strain.
The coefficients of the expansion only depend on the material of the deformed body and on
its physical state.

**Notation:**

#. :math:`O(f(x))` is the mathematical *Big O*. In this law the limit :math:`e \to 0` is
   assumed.

**Conditions:**

#. The deformations are elastic.
#. The deformations are small, i.e. :math:`e << 1`.
#. This law features the expansion up to the third power of strain, higher terms can be added
   if needed.
"""

from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

pressure = Symbol("pressure", units.pressure)
"""
Pressure (or tension) in the deformed body.

Symbol:
    :code:`P`
"""

young_modulus = Symbol("young_modulus", units.pressure)
"""
Young's modulus of the body's material.

Symbol:
    :code:`E`
"""

second_coefficient = Symbol("second_coefficient", units.pressure)
"""
Coefficient at the second power of strain in the expansion.

Symbol:
    :code:`A`
"""

third_coefficient = Symbol("third_coefficient", units.pressure)
"""
Coefficient at the third power of strain in the expansion.

Symbol:
    :code:`B`
"""

strain = Symbol("strain", dimensionless)
"""
Strain, or engineering normal strain to be exact, of the body.

Symbol:
    :code:`e`
"""

law = Eq(
    pressure,
    young_modulus * strain + second_coefficient * strain**2 + third_coefficient * strain**3,
)
r"""
:code:`P = E * e + A * e^2 + B * e^3 + O(e^4)`

Latex:
    .. math::
        P = E e + A e^2 + B e^3 + O(e^4)
"""


@validate_input(
    young_modulus_=young_modulus,
    second_coefficient_=second_coefficient,
    third_coefficient_=third_coefficient,
    strain_=strain,
)
@validate_output(pressure)
def calculate_pressure(
    young_modulus_: Quantity,
    second_coefficient_: Quantity,
    third_coefficient_: Quantity,
    strain_: float,
) -> Quantity:
    result = law.rhs.subs({
        young_modulus: young_modulus_,
        second_coefficient: second_coefficient_,
        third_coefficient: third_coefficient_,
        strain: strain_,
    })
    return Quantity(result)
