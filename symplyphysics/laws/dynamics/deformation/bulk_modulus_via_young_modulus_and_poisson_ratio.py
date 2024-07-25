r"""
Bulk modulus via Young modulus and Poisson ratio
================================================

Suppose a uniform isotropic body is subjected to bulk compression, i.e. forces are applied
to it from all its sides. Then the bulk modulus is the proportionality coefficient between
relative volume change of the body and the pressure inside of it. It is proportional to
the Young modulus of the material and also depends on its Poisson ratio.

**Conditions:**

#. :math:`\nu < \frac{1}{2}`, see the formula for elastic energy density.
..
    TODO add link to source file
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

bulk_modulus = Symbol("bulk_modulus", units.pressure)
"""
Bulk modulus of the material.

Symbol:
    :code:`K`
"""

young_modulus = Symbol("young_modulus", units.pressure)
"""
Young modulus of the material.

Symbol:
    :code:`E`
"""

poisson_ratio = Symbol("poisson_ratio", dimensionless)
r"""
Poisson ratio of the material.

Symbol:
    :code:`nu`

Latex:
    :math:`\nu`
"""

law = Eq(bulk_modulus, young_modulus / (3 * (1 - 2 * poisson_ratio)))
r"""
:code:`K = E / (3 * (1 - 2 * nu))`

Latex:
    .. math::
        K = \frac{E}{3 (1 - 2 \nu)}
"""


@validate_input(
    young_modulus_=young_modulus,
    poisson_ratio_=poisson_ratio,
)
@validate_output(bulk_modulus)
def calculate_bulk_modulus(
    young_modulus_: Quantity,
    poisson_ratio_: float,
) -> Quantity:
    result = law.rhs.subs({
        young_modulus: young_modulus_,
        poisson_ratio: poisson_ratio_,
    })
    return Quantity(result)
