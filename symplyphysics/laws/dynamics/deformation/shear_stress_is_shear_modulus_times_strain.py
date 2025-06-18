"""
Shear stress is shear modulus times strain
==========================================

When an object is under a shearing stress, the shear stress applied to the object is proportional to the shear modulus of the object and the angle formed by the face of the object perpendicular to the shearing force before (:math:`AB`) and after the deformation (:math:`AB'`), see **Note** for reference.

**Conditions:**

#. The deformation is elastic (reversible).

**Notes:**

#. For a visual representation of shear stress visit `this link <https://mechcontent.com/wp-content/uploads/2022/10/shear-deformation-in-object.webp>`__.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Shear_stress#Pure>`__.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.core.symbols.quantities import scale_factor

shear_stress = symbols.shear_stress
"""
:symbols:`shear_stress`
"""

shear_modulus = symbols.shear_modulus
"""
:symbols:`shear_modulus`.
"""

shear_strain = symbols.engineering_shear_strain
"""
:symbols:`engineering_shear_strain`.
"""

law = Eq(shear_stress, shear_modulus * shear_strain)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    shear_modulus_=shear_modulus,
    shear_strain_=shear_strain,
)
@validate_output(shear_stress)
def calculate_shear_stress(
    shear_modulus_: Quantity,
    shear_strain_: Quantity | float,
) -> Quantity:
    result = law.rhs.subs({
        shear_modulus: shear_modulus_,
        shear_strain: scale_factor(shear_strain_),
    })
    return Quantity(result)
