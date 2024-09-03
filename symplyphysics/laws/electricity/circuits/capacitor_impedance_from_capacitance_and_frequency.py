"""
Capacitor impedance from capacitance and frequency
==================================================

The impedance of the capacitor is a purely imaginary reactive impedance which
is inversely proportional to its capacitance. Note that the resistance of the
capacitor is infinite, i.e. it is considered to be an opened connection.
"""

from sympy import (I, Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, angle_type)

impedance = Symbol("impedance", units.impedance)
"""
Impedance of the capacitor.
"""

angular_frequency = Symbol("angular_frequency", angle_type / units.time)
r"""
Angular frequency of the alternating current.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

capacitance = Symbol("capacitance", units.capacitance)
"""
Capacitance of the capacitor.

Symbol:
    :code:`C`
"""

law = Eq(impedance, -I / (angular_frequency * capacitance))
r"""
:code:`Z = -i / (w * C)`

Latex:
    .. math::
        Z = -\frac{i}{\omega C}
"""


@validate_input(capacitance_=capacitance, circular_frequency_=angular_frequency)
@validate_output(impedance)
def calculate_impedance(capacitance_: Quantity, circular_frequency_: Quantity) -> Quantity:
    result_impedance_expr = solve(law, impedance, dict=True)[0][impedance]
    result_expr = result_impedance_expr.subs({
        capacitance: capacitance_,
        angular_frequency: circular_frequency_
    })
    return Quantity(result_expr)
