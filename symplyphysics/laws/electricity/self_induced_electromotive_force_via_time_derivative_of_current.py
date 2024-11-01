"""
Self-induced electromotive force via time derivative of current
===============================================================

Expression for the self-induced emf can be derived from the Faraday's law featuring
the time derivative of current flowing through the circuit.
"""

from sympy import Eq, Derivative
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
)
from symplyphysics.core.geometry.line import two_point_function, Point2D

time = symbols.time
"""
:symbols:`time`.
"""

electromotive_force = clone_as_function(symbols.electromotive_force, [time])
"""
Self-induced :symbols:`electromotive_force`.
"""

inductance = symbols.inductance
"""
:symbols:`inductance` of the circuit.
"""

current = clone_as_function(symbols.current, [time])
"""
:symbols:`current` in the circuit.
"""

law = Eq(electromotive_force(time), -1 * inductance * Derivative(current(time), time))
"""
:laws:symbol::

:laws:latex::
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
    applied_definition = law.subs({current(time): current_function_, inductance: inductance_})
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
