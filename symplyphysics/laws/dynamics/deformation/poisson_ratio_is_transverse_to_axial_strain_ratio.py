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
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

poisson_ratio = symbols.poisson_ratio
r"""
:symbols:`poisson_ratio` of the material of the deformed body.
"""

transverse_strain = clone_as_symbol(
    symbols.strain,
    display_symbol="e_transverse",
    display_latex="e_\\text{transverse}",
)
r"""
:symbols:`strain` in the transverse (lateral) direction relative to the deforming force.
"""

axial_strain = clone_as_symbol(
    symbols.strain,
    display_symbol="e_axial",
    display_latex="e_\\text{axial}",
)
r"""
:symbols:`strain` in the axial direction, i.e. parallel to the deforming force.
"""

law = Eq(poisson_ratio, -1 * transverse_strain / axial_strain)
"""
:laws:symbol::

:laws:latex::
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
