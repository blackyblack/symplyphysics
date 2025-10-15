"""
Mechanical work from force and displacement
===========================================

Work is measured result of force applied. Mechanical work is the only reason for the object energy
to be changed. Work is a scalar value equal to the dot product of force and displacement vectors.

**Conditions:**

#. Force must be constant during the body's displacement. This is possible when the displacement is
   infinitesimally small. For the case non-constant force and/or finite displacement, refer to
   :ref:`Mechanical work is line integral of force`.

**Notation:**

#. :math:`\\left( \\vec a, \\vec b \\right)` (:code:`dot(a, b)`) is the dot product between
   vectors :math:`\\vec a` and :math:`\\vec b`.

..
    TODO: Rename file to match heading
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorDot
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

work = symbols.work
"""
Mechanical :symbols:`work` done by the :attr:`~force` to displace the body.
"""

force = clone_as_vector_symbol(symbols.force)
"""
Vector of :symbols:`force` exerted on the body.
"""

displacement = clone_as_vector_symbol(symbols.distance)
"""
Vector denoting the change of the position vector of the body when it moves from one point in
space to another. Also see :symbols:`distance`.
"""

law = Eq(work, VectorDot(force, displacement))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(force_=force, distance_=displacement)
@validate_output(work)
def calculate_work(
    force_: QuantityCoordinateVector,
    distance_: QuantityCoordinateVector,
) -> Quantity:
    work_expr = law.rhs
    work_value = work_expr.subs({
        force: force_,
        displacement: distance_,
    })

    return Quantity(work_value)


# UNIQUE_LAW_ID: 241
