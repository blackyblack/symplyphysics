"""
Energy stored in inductor via inductance and current
====================================================

Inductors store energy in the magnetic field when the current is flowing
through it.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

energy = Symbol("energy", units.energy)
"""
Energy stored in the magnetic field.

Symbol:
    :code:`W`
"""

inductance = Symbol("inductance", units.inductance)
"""
Inductance of the inductor.

Symbol:
    :code:`L`
"""

current = Symbol("current", units.current)
"""
Current flowing through the inductor.

Symbol:
    :code:`I`
"""

law = Eq(energy, inductance * current**2 / 2)
r"""
:code:`W = 1/2 * L * I^2`

Latex:
    .. math::
        W = \frac{1}{2} L I^2
"""


@validate_input(inductance_=inductance, current_=current)
@validate_output(energy)
def calculate_accumulated_energy(inductance_: Quantity, current_: Quantity) -> Quantity:
    result_energy_expr = solve(law, energy, dict=True)[0][energy]
    result_expr = result_energy_expr.subs({inductance: inductance_, current: current_})
    return Quantity(result_expr)
