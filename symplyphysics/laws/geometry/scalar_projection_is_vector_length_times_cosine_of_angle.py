"""
Scalar projection is vector length times cosine of angle
========================================================

The **scalar projection** of a vector :math:`\\vec a` onto a vector :math:`\\vec b` can be found
if the vector_length of :math:`\\vec a` and the angle between :math:`\\vec a` and :math:`\\vec b` are
known.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Scalar_projection>`__.
"""

from sympy import Eq, solve, cos
from symplyphysics import Quantity, validate_input, symbols, Symbol
from symplyphysics.core.quantity_decorator import validate_output_same
from symplyphysics.core.symbols.quantities import scale_factor
from symplyphysics.core.dimensions import any_dimension

angle = symbols.angle
"""
:symbols:`angle` between :math:`\\vec a` and :math:`\\vec b`.
"""

vector_length = Symbol("a", any_dimension)
"""
Length of the projected vector :math:`\\vec a`.
"""

projection = Symbol("s", any_dimension)
"""
Length of the projection of :math:`\\vec a` onto :math:`\\vec b`.
"""

law = Eq(projection, vector_length * cos(angle))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(angle_=angle)
@validate_output_same("vector_length_")
def calculate_projection(vector_length_: Quantity, angle_: Quantity | float) -> Quantity:
    result_projection_expr = solve(law, projection, dict=True)[0][projection]
    #HACK: sympy angles are always in radians
    angle_radians = scale_factor(angle_)
    result_expr = result_projection_expr.subs({
        vector_length: vector_length_,
        angle: angle_radians,
    })
    return Quantity(result_expr)


# UNIQUE_LAW_ID: 180
