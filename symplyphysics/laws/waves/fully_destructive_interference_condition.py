"""
Fully destructive interference condition
========================================

The interference of two waves is said to be fully destructive when the amplitude of the
resulting wave is zero, i.e. the two waves cancel each other out. In that case, they must
be an odd number of half-cycles out of phase between each other.
"""

from sympy import Eq, pi
from symplyphysics import (
    angle_type,
    dimensionless,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
)

phase_shift = Symbol("phase_shift", angle_type, real=True)
r"""
Phase shift between interfering waves.

Symbol:
    :code:`phi`

Latex:
    :math:`\varphi`
"""

integer_factor = Symbol("integer_factor", dimensionless, integer=True)
"""
Integer factor.

Symbol:
    :code:`n`
"""

law = Eq(phase_shift, (1 + 2 * integer_factor) * pi)
r"""
:code:`phi = (1 + 2 * n) * pi`

Latex:
    .. math::
        \varphi = (1 + 2 n) \pi
"""


@validate_input(integer_factor_=integer_factor)
@validate_output(phase_shift)
def calculate_phase_shift(integer_factor_: int) -> Quantity:
    result = law.rhs.subs(integer_factor, integer_factor_)
    return Quantity(result)
