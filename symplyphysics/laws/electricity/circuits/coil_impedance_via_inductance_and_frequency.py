"""
Coil impedance via inductance and frequency
===========================================

The impedance of a coil depends on its inductance and the angular frequency of the alternating
current. It is a purely imaginary quantity since the resistance of a coil is zero.
"""

from sympy import (I, Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, angle_type)

impedance = Symbol("impedance", units.impedance)
"""
Impedance of the coil.

Symbol:
    :code:`Z`
"""

angular_frequency = Symbol("angular_frequency", angle_type / units.time)
r"""
Angular frequency of the current.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

inductance = Symbol("inductance", units.inductance)
"""
Coil inductance.

Symbol:
    :code:`L`
"""

law = Eq(impedance, I * angular_frequency * inductance)
r"""
:code:`Z = i * w * L`

Latex:
    .. math::
        Z = i \omega L
"""


@validate_input(inductance_=inductance, circular_frequency_=angular_frequency)
@validate_output(impedance)
def calculate_impedance(inductance_: Quantity, circular_frequency_: Quantity) -> Quantity:
    result_impedance_expr = solve(law, impedance, dict=True)[0][impedance]
    result_expr = result_impedance_expr.subs({
        inductance: inductance_,
        angular_frequency: circular_frequency_
    })
    return Quantity(result_expr)
