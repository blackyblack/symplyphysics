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
#. The deformations are small, i.e. :math:`e \ll 1`.
#. This law features the expansion up to the third power of strain, higher terms can be added
   if needed.
"""

from sympy import Eq, O
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

pressure = symbols.pressure
"""
:symbols:`pressure` (or tension) in the deformed body.
"""

young_modulus = symbols.young_modulus
"""
:symbols:`young_modulus` of the body's material.
"""

second_coefficient = clone_as_symbol(symbols.pressure, display_symbol="A", display_latex="A")
"""
Coefficient at the second power of strain in the expansion.
"""

third_coefficient = clone_as_symbol(symbols.pressure, display_symbol="B", display_latex="B")
"""
Coefficient at the third power of strain in the expansion.
"""

strain = symbols.engineering_normal_strain
"""
:symbols:`engineering_normal_strain` of the body.
"""

law = Eq(
    pressure,
    young_modulus * strain + second_coefficient * strain**2 + third_coefficient * strain**3 +
    O(strain**4),
)
"""
.. only:: comment

    Big O should be evaluated by SymPy

:laws:symbol::

:laws:latex::

:laws:sympy-eval::
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
    result = law.rhs.removeO().subs({
        young_modulus: young_modulus_,
        second_coefficient: second_coefficient_,
        third_coefficient: third_coefficient_,
        strain: strain_,
    })
    return Quantity(result)
