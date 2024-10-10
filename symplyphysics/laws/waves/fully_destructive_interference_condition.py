"""
Fully destructive interference condition
========================================

The interference of two waves is said to be fully destructive when the amplitude of the
resulting wave is zero, i.e. the two waves cancel each other out. In that case, they must
be an odd number of half-cycles out of phase between each other.
"""

from sympy import Eq, pi
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

phase_shift = symbols.phase_shift
"""
:symbols:`phase_shift` between interfering waves.
"""

integer_factor = symbols.whole_number
"""
Integer factor. See :symbols:`whole_number`.
"""

law = Eq(phase_shift, (1 + 2 * integer_factor) * pi)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(integer_factor_=integer_factor)
@validate_output(phase_shift)
def calculate_phase_shift(integer_factor_: int) -> Quantity:
    result = law.rhs.subs(integer_factor, integer_factor_)
    return Quantity(result)
