"""
Torque of twisting force
========================

Torque is a turning or twisting action on a body about a rotation axis due to a force. It is a
pseudovector defined as a cross product of the force vector and the position vector of the point
of force application.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Torque#Definition_and_relation_to_other_physical_quantities>`__.
"""

from sympy import Eq
from symplyphysics import (symbols, validate_input, validate_output)

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorCross
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector
from symplyphysics.core.experimental.solvers import solve_for_vector

torque = clone_as_vector_symbol(symbols.torque)
"""
Pseudovector of :symbols:`torque` due to twisting :attr:`~force`.
"""

force = clone_as_vector_symbol(symbols.force)
"""
Vector of :symbols:`force` exerted on the body.
"""

position_vector = clone_as_vector_symbol(symbols.distance_to_origin)
"""
Position vector of the body. See :symbols:`distance_to_origin`.
"""

law = Eq(torque, VectorCross(position_vector, force))
"""
:laws:symbol::

:laws:latex::
"""

# This law is actually the definition of torque in classical mechanics


@validate_input(force_=force, position_=position_vector)
@validate_output(torque)
def calculate_torque(
    position_: QuantityCoordinateVector,
    force_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = solve_for_vector(law, torque).subs({
        force: force_,
        position_vector: position_,
    })

    return QuantityCoordinateVector.from_expr(result)
