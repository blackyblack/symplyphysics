"""
Tensile stress is Young's modulus times strain
==============================================

When an object is under tension or compression, the stress is related to the strain via the
Young's modulus.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Young%27s_modulus#Definition>`__.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

stress = symbols.stress
"""
:symbols:`stress` on the object.
"""

young_modulus = symbols.young_modulus
"""
:symbols:`young_modulus` of the material of the object.
"""

engineering_normal_strain = symbols.engineering_normal_strain
"""
:symbols:`engineering_normal_strain` of the deformed body.
"""

law = Eq(stress, young_modulus * engineering_normal_strain)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    young_modulus_=young_modulus,
    engineering_normal_strain_=engineering_normal_strain,
)
@validate_output(stress)
def calculate_tensile_stress(
    young_modulus_: Quantity,
    engineering_normal_strain_: float,
) -> Quantity:
    result = law.rhs.subs({
        young_modulus: young_modulus_,
        engineering_normal_strain: engineering_normal_strain_,
    })
    return Quantity(result)
