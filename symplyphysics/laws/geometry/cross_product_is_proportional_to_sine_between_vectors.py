"""
Cross product is proportional to sine of angle between vectors
==============================================================

The cross product of two vectors is a binary operation which produces a vector whose length is
proportional to the lengths of the given vectors and the sine of the angle between them.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Cross_product#Definition>`__.
"""

from sympy import Eq, sin, solve
from symplyphysics import Quantity, Symbol, symbols
from symplyphysics.core.dimensions import any_dimension
from symplyphysics.core.symbols.quantities import scale_factor

cross_product_length = Symbol(
    "norm(cross(u, v))",
    any_dimension,
    display_latex="\\left \\Vert \\left[ \\vec u, \\vec v \\right] \\right \\Vert",
    positive=True,
)
"""
Length of the cross product between :math:`\\vec u` and :math:`\\vec v`.
"""

first_vector_length = Symbol("u", any_dimension, positive=True)
"""
Length of :math:`\\vec u`.
"""

second_vector_length = Symbol("v", positive=True)
"""
Length of :math:`\\vec v`.
"""

angle_between_vectors = symbols.angle
"""
:symbols:`angle` between :math:`\\vec u` and :math:`\\vec v`.
"""

law = Eq(
    cross_product_length,
    first_vector_length * second_vector_length * sin(angle_between_vectors),
)
"""
:laws:symbol::

:laws:latex::
"""


def calculate_sine_between_vectors(
    cross_product_norm_: Quantity | float,
    vector_left_norm_: Quantity | float,
    vector_right_norm_: Quantity | float,
) -> Quantity:
    cross_product_norm_ = scale_factor(cross_product_norm_)
    vector_left_norm_ = scale_factor(vector_left_norm_)
    vector_right_norm_ = scale_factor(vector_right_norm_)
    result = solve(law, sin(angle_between_vectors))[0].subs({
        cross_product_length: cross_product_norm_,
        first_vector_length: vector_left_norm_,
        second_vector_length: vector_right_norm_,
    })
    return Quantity(result)


# UNIQUE_LAW_ID: 181
