"""
Solution to the exponential decay equation
==========================================

The solution to the exponential decay equation is the product of the initial quantity
and the the ratio of the current time to the half-life of the quantity, raised to the
power of 2. In other words, for every half-life that passes, the quantity decays by a
factor of 2.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Exponential_decay#>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

final_quantity = symbols.any_quantity
"""
Quantity that still remains and has not decayed after :attr:`~time` :math:`t`.
"""

initial_quantity = clone_as_symbol(symbols.any_quantity, subscript="0")
"""
Initial quantity that will decay.
"""

half_life = symbols.half_life
"""
:symbols:`half_life` of the decaying quantity.
"""

time = symbols.time
"""
:symbols:`time`.
"""

law = Eq(final_quantity, initial_quantity * 2**(-1 * time / half_life))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(half_life_=half_life,
    decay_time_=time)
def calculate_number_of_cores(number_of_cores_initial_: int, half_life_: Quantity,
    decay_time_: Quantity) -> int:
    if number_of_cores_initial_ < 0:
        raise ValueError("Number of cores cannot be negative")
    result_expr = solve(law, final_quantity, dict=True)[0][final_quantity]
    result_expr = result_expr.subs({
        initial_quantity: number_of_cores_initial_,
        half_life: half_life_,
        time: decay_time_
    })
    return int(convert_to_float(result_expr))
