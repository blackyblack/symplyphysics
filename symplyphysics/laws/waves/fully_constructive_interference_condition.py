"""
Fully constructive interference condition
=========================================

The interference of two waves is said to be *fully constructive* when the amplitude of the
resulting wave is precisely the sum of the amplitudes of the comprising waves. In that case,
there must be an even number of half-cycles out of phase between each other.
"""

from sympy import Eq, pi
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

phase_shift = symbols.phase_shift
r"""
:symbols:`phase_shift` between interfering waves.
"""

integer_factor = symbols.whole_number
"""
Integer factor. See :symbols:`whole_number`
"""

law = Eq(phase_shift, 2 * pi * integer_factor)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(integer_factor_=integer_factor)
@validate_output(phase_shift)
def calculate_phase_shift(integer_factor_: int) -> Quantity:
    result = law.rhs.subs(integer_factor, integer_factor_)
    return Quantity(result)
