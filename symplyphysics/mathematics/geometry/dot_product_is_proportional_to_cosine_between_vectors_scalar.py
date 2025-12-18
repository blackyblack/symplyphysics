"""
Dot product is proportional to cosine of angle between vectors
==============================================================

The dot product of two vectors is a scalar binary operation that can be defined as the product of
the norms of the vectors and the cosine of the angle between them. Also see the :ref:`vector law
<Dot product is proportional to cosine of angle between vectors (vector)>`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Dot_product#Geometric_definition>`__.
"""

from sympy import Eq, cos, solve
from symplyphysics import Quantity, Symbol, symbols
from symplyphysics.core.dimensions import any_dimension
from symplyphysics.core.symbols.quantities import scale_factor

dot_product = Symbol(
    "dot(u, v)",
    any_dimension,
    display_latex="\\left( \\vec u, \\vec v \\right)",
    real=True,
)
"""
Dot product between :math:`\\vec u` and :math:`\\vec v`.
"""

first_vector_length = Symbol("u", any_dimension, positive=True)
"""
Length of :math:`\\vec u`.
"""

second_vector_length = Symbol("v", any_dimension, positive=True)
"""
Length of :math:`\\vec v`.
"""

angle_between_vectors = symbols.angle
"""
:symbols:`angle` between :math:`\\vec u` and :math:`\\vec v`.
"""

law = Eq(dot_product, first_vector_length * second_vector_length * cos(angle_between_vectors))
"""
:laws:symbol::

:laws:latex::
"""


def calculate_cosine_between_vectors(
    dot_product_: Quantity | float,
    vector_left_norm_: Quantity | float,
    vector_right_norm_: Quantity | float,
) -> Quantity:
    dot_product_ = scale_factor(dot_product_)
    vector_left_norm_ = scale_factor(vector_left_norm_)
    vector_right_norm_ = scale_factor(vector_right_norm_)
    result = solve(law, cos(angle_between_vectors))[0].subs({
        dot_product: dot_product_,
        first_vector_length: vector_left_norm_,
        second_vector_length: vector_right_norm_,
    })
    return Quantity(result)
