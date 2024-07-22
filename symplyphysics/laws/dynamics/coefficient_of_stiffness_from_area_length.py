"""
Coefficient of stiffness from area and length
=========================================

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

coefficient_of_stiffness = Symbol("coefficient_of_stiffness", units.force / units.length)
"""
The stiffness coefficient of the material.

Symbol:
    k
"""

young_modulus = Symbol("young_modulus", units.pressure)
"""
The Young's modulus of the material.

Symbol:
    E
"""

area = Symbol("area", units.area)
"""
The cross-sectional area of the object.

Symbol:
    A
"""

length = Symbol("length", units.length)
"""
The length of the object.

Symbol:
    l
"""

law = Eq(coefficient_of_stiffness, young_modulus * area / length)
r"""
k = E * A / l

Latex:
    :math:`k = E \frac{A}{l}`
"""


@validate_input(module_of_young_=young_modulus, area_=area, length_=length)
@validate_output(coefficient_of_stiffness)
def calculate_coefficient_of_stiffness(module_of_young_: Quantity, area_: Quantity,
    length_: Quantity) -> Quantity:
    result_expr = solve(law, coefficient_of_stiffness, dict=True)[0][coefficient_of_stiffness]
    result_expr = result_expr.subs({
        young_modulus: module_of_young_,
        area: area_,
        length: length_
    })
    return Quantity(result_expr)
