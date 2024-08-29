"""
Self-induced emf via time derivative of current
===============================================

Expression for the self-induced emf can be derived from the Faraday's law featuring
the time derivative of current flowing through the circuit.
"""

from sympy import (Eq, Derivative)
from symplyphysics import (units, Quantity, Function, Symbol, validate_input,
    validate_output)
from symplyphysics.core.geometry.line import two_point_function, Point2D

electromotive_force = Function("electromotive_force", units.voltage)
r"""
Self-induced electromotive force.

Symbol:
    :code:`E(t)`

Latex:
    :math:`\mathcal{E}(t)`
"""

inductance = Symbol("inductance", units.inductance)
"""
Inductance of the circuit.

Symbol:
    :code:`L`
"""

current = Function("current", units.current)
"""
Current in the circuit.

Symbol:
    :code:`I(t)`
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

law = Eq(electromotive_force(time), -1 * inductance * Derivative(current(time), time))
r"""
:code:`E(t) = -1 * L * Derivative(I(t), t)`

Latex:
    .. math::
        \mathcal{E}(t) = -L \frac{d I}{d t}
"""


@validate_input(inductance_=inductance, current_start_=current, current_end_=current, time_=time)
@validate_output(electromotive_force)
def calculate_voltage(inductance_: Quantity, current_start_: Quantity, current_end_: Quantity,
    time_: Quantity) -> Quantity:
    current_function_ = two_point_function(
        Point2D(0, current_start_),
        Point2D(time_, current_end_),
        time,
    )
    applied_definition = law.subs({
        current(time): current_function_,
        inductance: inductance_
    })
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
