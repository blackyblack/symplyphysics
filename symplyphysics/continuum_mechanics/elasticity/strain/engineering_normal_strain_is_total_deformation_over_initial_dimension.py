"""
Engineering normal strain is total deformation over initial dimension
=====================================================================

*Engineering normal strain*, also called *Cauchy strain*, is expressed as the ratio of total
deformation to the initial dimension of a material body on which forces are applied.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Strain_(mechanics)#Engineering_strain>`__.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
)

engineering_normal_strain = symbols.engineering_normal_strain
"""
:symbols:`engineering_normal_strain` of the body.
"""

total_deformation = symbols.deformation
"""
Total :symbols:`deformation` of the body.
"""

initial_dimension = symbols.length
"""
Initial dimension of the body, e.g. initial :symbols:`length`.
"""

law = Eq(engineering_normal_strain, total_deformation / initial_dimension)
"""
:laws:symbol::

:laws:latex::
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
