"""
Engineering normal strain is total deformation over initial dimension
=====================================================================

*Engineering normal strain*, also called *Cauchy strain*, is expressed as the ratio of total
deformation to the initial dimension of a material body on which forces are applied.
"""

from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    convert_to_float,
)

engineering_normal_strain = Symbol("engineering_normal_strain", dimensionless)
r"""
*Engineering normal strain*, also known as *engineering extensional strain*, is equal to
the *relative elongation* of the body.

Symbol:
    :code:`e`
"""

total_deformation = Symbol("total_deformation", units.length)
r"""
Total deformation of the body, which is the change in length in case of the engineering
normal strain.

Symbol:
    :code:`dl`

Latex:
    :math:`\Delta l`
"""

initial_dimension = Symbol("initial_dimension", units.length)
"""
Initial dimension of the body, e.g. initial length.

Symbol:
    :code:`l`
"""

law = Eq(engineering_normal_strain, total_deformation / initial_dimension)
r"""
:code:`e = dl / l`

Latex:
    .. math::
        e = \frac{\Delta l}{l}
"""


@validate_input(
    total_deformation_=total_deformation,
    initial_dimension_=initial_dimension,
)
@validate_output(engineering_normal_strain)
def calculate_engineering_normal_strain(
    total_deformation_: Quantity,
    initial_dimension_: Quantity,
) -> float:
    result = law.rhs.subs({
        total_deformation: total_deformation_,
        initial_dimension: initial_dimension_,
    })
    return convert_to_float(result)
