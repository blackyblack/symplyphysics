"""
Restoring torque due to twist of torsion pendulum
=================================================

A torsion pendulum is a simple harmonic oscillator consisting of a disk suspended by a wire.
Rotating the disk through an angle in either direction introduces a restoring torque.

**Links:**

#. `Wikipedia, scalar equation <https://en.wikipedia.org/wiki/Torsion_spring#Torsion_coefficient>`__.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols

from symplyphysics.core.experimental.solvers import solve_for_vector
from symplyphysics.core.experimental.vectors import clone_as_vector_symbol
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector

torque = clone_as_vector_symbol(symbols.torque)
"""
Pseudovector of :symbols:`torque` of the pendulum.
"""

torsion_stiffness = symbols.torsion_stiffness
"""
:symbols:`torsion_stiffness` of the pendulum.
"""

angular_displacement = clone_as_vector_symbol(symbols.angular_distance)
"""
Pseudovector of angular displacement. See :symbols:`angular_distance`.
"""

law = Eq(torque, -1 * torsion_stiffness * angular_displacement)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(torsion_constant_=torsion_stiffness, rotation_vector_=angular_displacement)
@validate_output(torque)
def calculate_torque(
    torsion_constant_: Quantity,
    rotation_vector_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    expr = solve_for_vector(law, torque)
    value = expr.subs({
        torsion_stiffness: torsion_constant_,
        angular_displacement: rotation_vector_,
    })

    return QuantityCoordinateVector.from_expr(value)


@validate_input(torsion_constant_=torsion_stiffness, torque_=torque)
@validate_output(angular_displacement)
def calculate_rotation_vector(
    torsion_constant_: Quantity,
    torque_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    expr = solve_for_vector(law, angular_displacement)
    value = expr.subs({
        torsion_stiffness: torsion_constant_,
        torque: torque_,
    })

    return QuantityCoordinateVector.from_expr(value)
