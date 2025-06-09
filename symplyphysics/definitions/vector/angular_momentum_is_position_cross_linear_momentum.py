r"""
Angular momentum is position cross linear momentum
==================================================

The pseudovector of *angular momentum* of a particle is the cross product of its
position vector and linear momentum. Unlike linear momentum, angular momentum depends
on where the origin of the coordinate system is chosen since it depends on the position
vector of the particle defined relative to that origin.

**Notation:**

#. :math:`\vec a \times \vec b` (:code:`cross(a, b)`) denotes a cross product between
   vectors :math:`\vec a` and :math:`\vec b`.

**Links:**

#. `Wikipedia, see second paragraph <https://en.wikipedia.org/wiki/Angular_momentum#>`__.
"""

from sympy import Eq, Expr
from symplyphysics import validate_input, validate_output, symbols

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorCross
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

position_vector = clone_as_vector_symbol(symbols.distance_to_origin)
"""
Displacement of the body relative to the origin of the reference frame. See
:symbols:`distance_to_origin`.
"""

linear_momentum = clone_as_vector_symbol(symbols.momentum)
"""
Vector of linear :symbols:`momentum`.
"""

angular_momentum = clone_as_vector_symbol(symbols.angular_momentum)
"""
Pseudovector of :symbols:`angular_momentum`.
"""

law = Eq(angular_momentum, VectorCross(position_vector, linear_momentum))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    position_vector_=position_vector,
    linear_momentum_=linear_momentum,
)
@validate_output(angular_momentum)
def calculate_angular_momentum(
    position_vector_: QuantityCoordinateVector,
    linear_momentum_: QuantityCoordinateVector,
) -> Expr:
    angular_momentum_value = law.rhs.subs({
        position_vector: position_vector_,
        linear_momentum: linear_momentum_,
    })

    return QuantityCoordinateVector.from_expr(angular_momentum_value)
