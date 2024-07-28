"""
Poisson ratio is transverse to axial strain ratio
=================================================

Poisson ratio is a physical quantity equal to the ratio of transverse strain
to axial (or longitudinal) strain in the direction of a deforming force.

**Notes:**

#. The sign convention is as follows: positive strain indicates extension and
   negative strain indicates contraction.
#. The Poisson ratio and the Young modulus are all that is needed to completely
   describe elastic properties of an isotropic material.
"""

from sympy import Eq
from symplyphysics import (
    dimensionless,
    Symbol,
    validate_input,
    validate_output,
    convert_to_float,
)

poisson_ratio = Symbol("poisson_ratio", dimensionless)
r"""
Poisson ratio of the material of the deformed body.

Symbol:
    :code:`nu`

Latex:
    :math:`\nu`
"""

transverse_strain = Symbol("transverse_strain", dimensionless)
r"""
Strain in the transverse (lateral) direction relative to the deforming force.

Symbol:
    :code:`e_transverse`

Latex:
    :math:`e_\text{transverse}`
"""

axial_strain = Symbol("axial_strain", dimensionless)
r"""
Strain in the axial direction, i.e. parallel to the deforming force.

Symbol:
    :code:`e_axial`

Latex:
    :math:`e_\text{axial}`
"""

law = Eq(poisson_ratio, -1 * transverse_strain / axial_strain)
r"""
:code:`nu = -1 * e_transverse / e_axial`

Latex:
    .. math::
        \nu = - \frac{e_\text{transverse}}{e_\text{axial}}
"""


@validate_input(
    transverse_strain_=transverse_strain,
    axial_strain_=axial_strain,
)
@validate_output(poisson_ratio)
def calculate_poisson_ratio(
    transverse_strain_: float,
    axial_strain_: float,
) -> float:
    result = law.rhs.subs({
        transverse_strain: transverse_strain_,
        axial_strain: axial_strain_,
    })
    return convert_to_float(result)
