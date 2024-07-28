"""
Tensile stress is Young's modulus times strain
==============================================

When an object is under tension or compression, the stress is related to the strain via the
Young's modulus.
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

stress = Symbol("stress", units.pressure)
r"""
Stress on the object.

Symbol:
    :code:`sigma`

Latex:
    :math:`\sigma`
"""

youngs_modulus = Symbol("youngs_modulus", units.pressure)
"""
Young's modulus of the material of the object.

Symbol:
    :code:`E`
"""

engineering_normal_strain = Symbol("engineering_normal_strain", dimensionless)
"""
Engineering normal strain of the deformed body.

Symbol:
    :code:`e`
"""

law = Eq(stress, youngs_modulus * engineering_normal_strain)
r"""
:code:`sigma = E * e`

Latex:
    .. math::
        \sigma = E e
"""


@validate_input(
    youngs_modulus_=youngs_modulus,
    engineering_normal_strain_=engineering_normal_strain,
)
@validate_output(stress)
def calculate_tensile_stress(
    youngs_modulus_: Quantity,
    engineering_normal_strain_: float,
) -> Quantity:
    result = law.rhs.subs({
        youngs_modulus: youngs_modulus_,
        engineering_normal_strain: engineering_normal_strain_,
    })
    return Quantity(result)
