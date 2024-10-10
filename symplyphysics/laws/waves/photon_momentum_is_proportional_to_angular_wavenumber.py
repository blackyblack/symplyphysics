"""
Photon momentum is proportional to angular wavenumber
=====================================================

The momentum of a photon is proportional to its wavenumber.

**Notation:**

#. :quantity_notation:`hbar`.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    quantities,
    symbols,
)

momentum = symbols.momentum
"""
:symbols:`momentum` of a photon.
"""

angular_wavenumber = symbols.angular_wavenumber
"""
:symbols:`angular_wavenumber` of a photon.
"""

law = Eq(momentum, quantities.hbar * angular_wavenumber)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(wavenumber_=angular_wavenumber)
@validate_output(momentum)
def calculate_momentum(wavenumber_: Quantity) -> Quantity:
    result_momentum_expr = solve(law, momentum, dict=True)[0][momentum]
    result_expr = result_momentum_expr.subs({angular_wavenumber: wavenumber_})
    return Quantity(result_expr)
