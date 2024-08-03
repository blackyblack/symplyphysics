"""
Coefficient of stiffness from area and length
=============================================

The *stiffness* coefficient depends on the Young's modulus, cross-sectional area and length.
Young's modulus is a tabular value that is different for each material.
"""

from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

stiffness = Symbol("stiffness", units.force / units.length)
"""
The stiffness coefficient of the material.

Symbol:
    :code:`k`
"""

young_modulus = Symbol("young_modulus", units.pressure)
"""
The Young's modulus of the material.

Symbol:
    :code:`E`
"""

area = Symbol("area", units.area)
"""
The cross-sectional area of the object.

Symbol:
    :code:`A`
"""

length = Symbol("length", units.length)
"""
The length of the object.

Symbol:
    :code:`l`
"""

law = Eq(stiffness, young_modulus * area / length)
r"""
:code:`k = E * A / l`

Latex:
    :math:`k = E \frac{A}{l}`
"""


@validate_input(module_of_young_=young_modulus, area_=area, length_=length)
@validate_output(stiffness)
def calculate_coefficient_of_stiffness(module_of_young_: Quantity, area_: Quantity,
    length_: Quantity) -> Quantity:
    result_expr = solve(law, stiffness, dict=True)[0][stiffness]
    result_expr = result_expr.subs({young_modulus: module_of_young_, area: area_, length: length_})
    return Quantity(result_expr)
