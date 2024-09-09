"""
Time constant of resistor-capacitor circuit
===========================================

The time constant of an RC circuit is the time it takes for the current to become
:math:`e` times smaller its original value. Also see `exponential decay
<https://en.wikipedia.org/wiki/Exponential_decay>`_.
"""

from sympy import Eq
from symplyphysics import units, Quantity, Symbol, validate_input, validate_output

time_constant = Symbol("time_constant", units.time)
r"""
Time constant of the circuit.

Symbol:
    :code:`tau`

Latex:
    :math:`\tau`
"""

resistance = Symbol("resistance", units.impedance)
"""
Resistance of the resistor.

Symbol:
    :code:`R`
"""

capacitance = Symbol("capacitance", units.capacitance)
"""
Capacitance of the capacitor.

Symbol:
    :code:`C`
"""

law = Eq(time_constant, resistance * capacitance)
r"""
:code:`tau = R * C`

Latex:
    .. math::
        \tau = R C
"""


@validate_input(resistance_=resistance, capacitance_=capacitance)
@validate_output(time_constant)
def calculate_time_constant(resistance_: Quantity, capacitance_: Quantity) -> Quantity:
    result = law.rhs.subs({
        resistance: resistance_,
        capacitance: capacitance_,
    })
    return Quantity(result)
