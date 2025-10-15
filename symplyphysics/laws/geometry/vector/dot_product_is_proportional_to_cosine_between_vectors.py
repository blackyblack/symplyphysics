"""
Dot product is proportional to cosine of angle between vectors (vector)
=======================================================================

The dot product of two vectors is a scalar binary operation that can be defined as the product of
the norms of the vectors and the cosine of the angle between them. Also see the :ref:`scalar law
<Dot product is proportional to cosine of angle between vectors>`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Dot_product#Geometric_definition>`__.
"""

from sympy import solve, Eq, cos

from symplyphysics import Quantity, symbols
from symplyphysics.core.dimensions import any_dimension

from symplyphysics.core.experimental.vectors import VectorSymbol, VectorDot, VectorNorm
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

first_vector = VectorSymbol("u", any_dimension)
"""
First vector.
"""

second_vector = VectorSymbol("v", any_dimension)
"""
Second vector.
"""

angle_between_vectors = symbols.angle
"""
:symbols:`angle` between :attr:`~first_vector` and :attr:`~second_vector`.
"""

law = Eq(
    VectorDot(first_vector, second_vector),
    VectorNorm(first_vector) * VectorNorm(second_vector) * cos(angle_between_vectors),
)
"""
:laws:symbol::

:laws:latex::
"""


def calculate_cosine_between_vectors(
    vector_left_: QuantityCoordinateVector,
    vector_right_: QuantityCoordinateVector,
) -> Quantity:
    result = solve(law, cos(angle_between_vectors))[0].subs({
        first_vector: vector_left_,
        second_vector: vector_right_,
    }).doit()

    return Quantity(result)


# UNIQUE_LAW_ID: 184
