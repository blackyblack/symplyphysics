"""
Fully constructive interference condition
=========================================

The interference of two waves is said to be *fully constructive* when the amplitude of the
resulting wave is precisely the sum of the amplitudes of the comprising waves.
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

law = Eq(phase_shift, 2 * pi * integer_factor)
r"""
:code:`phi = 2 * pi * n`

Latex:
    .. math::
        \varphi = 2 \pi n
"""

@validate_input(integer_factor_=integer_factor)
@validate_output(phase_shift)
def calculate_phase_shift(integer_factor_: int) -> Quantity:
    result = law.rhs.subs(integer_factor, integer_factor_)
    return Quantity(result)
