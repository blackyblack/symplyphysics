"""
Initial momentum equals final momentum
======================================

If there is no external force applied to system of objects, the summary momentum of this
system remains constant during and after any interactions between objects. See
:ref:`Momentum is constant`.

**Conditions:**

#. The system is closed.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Momentum#Conservation>`__.

..
    TODO: rename file
"""

from sympy import Eq, solve, dsolve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.conservation import momentum_of_colliding_objects_is_constant as constant_momentum

initial_time = clone_as_symbol(symbols.time, subscript="0")
"""
Initial :symbols:`time`.
"""

final_time = clone_as_symbol(symbols.time, subscript="1")
"""
Final :symbols:`time`.
"""

momentum = clone_as_function(symbols.momentum, [symbols.time])
"""
:symbols:`momentum` as a function of :symbols:`time`.
"""

law = Eq(momentum(final_time), momentum(initial_time))
"""
:laws:symbol::

:laws:latex::
"""

# Derive the same law from constant momentum

## dsolve() shows that solution is constant C1
_dsolved = dsolve(constant_momentum.law, constant_momentum.momentum(constant_momentum.time))

_energy_before_eq = _dsolved.subs(constant_momentum.time, initial_time)
_energy_before_eq = _energy_before_eq.subs(constant_momentum.momentum(initial_time),
    momentum(initial_time))
_energy_after_eq = _dsolved.subs(constant_momentum.time, final_time)
_energy_after_eq = _energy_after_eq.subs(constant_momentum.momentum(final_time),
    momentum(final_time))

## Show that when energy is constant, energy_before equals to energy_after
_energy_after_solved = solve([_energy_after_eq, _energy_before_eq], (momentum(final_time), "C1"),
    dict=True)[0][momentum(final_time)]
assert expr_equals(_energy_after_solved, law.rhs)


@validate_input(momentum_before_=momentum)
@validate_output(momentum)
def calculate_momentum_after(momentum_before_: Quantity) -> Quantity:
    solved = solve(law, momentum(final_time), dict=True)[0][momentum(final_time)]
    result_expr = solved.subs(momentum(initial_time), momentum_before_)
    return Quantity(result_expr)
