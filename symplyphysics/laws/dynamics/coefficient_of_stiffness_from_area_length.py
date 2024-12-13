"""
Coefficient of stiffness from area and length
=============================================

The *stiffness* coefficient depends on the Young's modulus, cross-sectional area and length.
Young's modulus is a tabular value that is different for each material.

**Links:**

#. `Wikipedia, first formula <https://en.wikipedia.org/wiki/Stiffness#Relationship_to_elasticity>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

stiffness = symbols.stiffness
"""
The :symbols:`stiffness` coefficient of the material.
"""

young_modulus = symbols.young_modulus
"""
The :symbols:`young_modulus` of the material.
"""

area = symbols.area
"""
The cross-sectional :symbols:`area` of the object.
"""

length = symbols.length
"""
The :symbols:`length` of the object.
"""

law = Eq(stiffness, young_modulus * area / length)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(module_of_young_=young_modulus, area_=area, length_=length)
@validate_output(stiffness)
def calculate_coefficient_of_stiffness(module_of_young_: Quantity, area_: Quantity,
    length_: Quantity) -> Quantity:
    result_expr = solve(law, stiffness, dict=True)[0][stiffness]
    result_expr = result_expr.subs({young_modulus: module_of_young_, area: area_, length: length_})
    return Quantity(result_expr)
