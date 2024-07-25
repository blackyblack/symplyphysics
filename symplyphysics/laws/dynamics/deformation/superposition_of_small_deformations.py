"""
Superposition of small deformations
===================================

If the deformations are small, the elastic properties of bodies stay the same under deformations.
If several forces act on a body, its total deformation can be computed from the deformations of
the individual forces, however these deformations must be small.

**Conditions:**

#. The deformations are small.
"""

from sympy import Eq
from symplyphysics import (
    dimensionless,
    Symbol,
    validate_input,
    validate_output,
    convert_to_float,
)

total_strain = Symbol("total_strain", dimensionless)
"""
Total strain of the body.

Symbol:
    :code:`e`
"""

first_strain = Symbol("first_strain", dimensionless)
r"""
Strain caused by one force.

Symbol:
    :code:`e_1`

Latex:
    :math:`e_1`
"""

second_strain = Symbol("second_strain", dimensionless)
r"""
Strain caused by another force.

Symbol:
    :code:`e_2`

Latex:
    :math:`e_2`
"""

law = Eq(total_strain, first_strain + second_strain)
r"""
:code:`e = e_1 + e_2`

Latex:
    .. math::
        e = e_1 + e_2
"""


@validate_input(
    first_strain_=first_strain,
    second_strain_=second_strain,
)
@validate_output(total_strain)
def calculate_total_strain(
    first_strain_: float,
    second_strain_: float,
) -> float:
    result = law.rhs.subs({
        first_strain: first_strain_,
        second_strain: second_strain_,
    })
    return convert_to_float(result)
