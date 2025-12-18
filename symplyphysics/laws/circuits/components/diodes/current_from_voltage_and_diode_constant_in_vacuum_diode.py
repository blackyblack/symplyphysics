"""
Current from voltage and diode constant in vacuum diode
=======================================================

The current-voltage characteristic of a vacuum diode is described by the
:math:`3/2`-power law (referred to as **Child's law** or **Child—Langmuir law**). The
diode constant in this law depends only on the relative position, shape and size of the
electrodes of the vacuum diode.

**Conditions:**

#. Electrons travel ballistically between electrodes, i.e. without scattering.
#. In the interelectrode region, the space charge of any ions is negligible.
#. The electrons have zero velocity at the cathode surface.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Space_charge#In_vacuum_(Child's_law)>`__.
"""

from sympy import Eq, solve, Rational
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

anode_current = symbols.current
"""
Anode :symbols:`current`.
"""

diode_constant = symbols.diode_constant
"""
:symbols:`diode_constant`.
"""

anode_voltage = clone_as_symbol(symbols.voltage, display_symbol="U_a", display_latex="U_\\text{a}")
"""
:symbols:`voltage` between cathode and anode.
"""

law = Eq(anode_current, diode_constant * anode_voltage**Rational(3, 2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(diode_constant_=diode_constant, voltage_=anode_voltage)
@validate_output(anode_current)
def calculate_current(diode_constant_: Quantity, voltage_: Quantity) -> Quantity:
    result_expr = solve(law, anode_current, dict=True)[0][anode_current]
    result_expr = result_expr.subs({
        diode_constant: diode_constant_,
        anode_voltage: voltage_,
    })
    return Quantity(result_expr)
