r"""
Photon momentum is proportional to angular wavenumber
=====================================================

The momentum of a photon is proportional to its wavenumber.

**Notation:**

#. :quantity_notation:`hbar`.
"""

from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    quantities,
)

momentum = Symbol("momentum", units.momentum)
"""
Momentum of a photon.

Symbol:
    :code:`p`
"""

angular_wavenumber = Symbol("angular_wavenumber", 1 / units.length)
"""
Angular wavenumber of a photon.

Symbol: 
    :code:`k`
"""

law = Eq(momentum, quantities.hbar * angular_wavenumber)
r"""
:code:`p = hbar * k`

Latex:
    .. math::
        p = \hbar k
"""


@validate_input(wavenumber_=angular_wavenumber)
@validate_output(momentum)
def calculate_momentum(wavenumber_: Quantity) -> Quantity:
    result_momentum_expr = solve(law, momentum, dict=True)[0][momentum]
    result_expr = result_momentum_expr.subs({angular_wavenumber: wavenumber_})
    return Quantity(result_expr)
