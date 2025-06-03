"""
Force is derivative of momentum (vector)
========================================

Newton's second law of motion can be generalized in terms of linear momentum. Precisely,
the net force exerted on a body is equal to the time derivative of the body's momentum.

**Notes:**

#. Works in relativistic mechanics as well as in classical mechanics.

#. See :ref:`scalar counterpart <Force is derivative of momentum>` of this law.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Momentum#Relation_to_force>`__.
"""

from sympy import Eq, Expr

from symplyphysics import units, Quantity, symbols
from symplyphysics.core.dimensions import assert_equivalent_dimension

from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector
from symplyphysics.core.experimental.vectors import VectorFunction, VectorDerivative
from symplyphysics.core.experimental.solvers import solve_for_vector

time = symbols.time
"""
:symbols:`time`.
"""

momentum = VectorFunction("p", (time,), dimension=units.momentum, display_latex="\\mathbf{p}")
"""
The magnitude of the :symbols:`momentum` of the body as a function of :attr:`~time`.
"""

force = VectorFunction("F", (time,), dimension=units.force, display_latex="\\mathbf{F}")
"""
The magnitude of the net :symbols:`force` exerted on the body as a function of :attr:`~time`.
"""

force_law = Eq(force(time), VectorDerivative(momentum(time), time))
"""
..
    Auto-printing of vector expressions has not been implemented yet.

:code:`F(t) = Derivative(p(t), t)`

Latex:
    .. math::
        \\mathbf{F}(t) = \\frac{d}{d t} \\mathbf{p}(t)
"""


def calculate_force(
    momentum_before_: QuantityCoordinateVector,
    momentum_after_: QuantityCoordinateVector,
    time_: Quantity,
) -> Expr:
    # TODO: fix `validate_input` to account for vector expressions

    assert_equivalent_dimension(momentum_before_.dimension, "momentum_before_", "calculate_force",
        momentum.dimension)
    assert_equivalent_dimension(momentum_after_.dimension, "momentum_after_", "calculate_force",
        momentum.dimension)
    assert_equivalent_dimension(time_, "time_", "calculate_force", time.dimension)

    momentum_ = (time / time_) * (momentum_after_ - momentum_before_)

    solved = solve_for_vector(force_law, force(time)).rhs
    force_ = solved.subs(momentum(time), momentum_).doit()
    force_ = QuantityCoordinateVector.combine(force_)

    assert_equivalent_dimension(force_.dimension, "result", "calculate_force", force.dimension)

    return force_
