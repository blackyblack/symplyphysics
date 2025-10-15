"""
Current is voltage over impedance
==================================

Alternating current flowing through an electric component is proportional to voltage applied to it
and inversely proportional to its impedance. This is the complex generalization of the :ref:`Ohm's
law <Current is voltage over resistance>`.

**Notes:**

#. This is a phenomenological equation.

**Conditions:**

#. The material must be ohmic, i.e. current must be proportional to voltage.

#. The AC circuit must be time-invariant.

**Links:**

#. `Wikipedia — Ohm's law <https://en.wikipedia.org/wiki/Ohm%27s_law#Reactive_circuits_with_time-varying_signals>`__.

#. `Univeristy of Maryland — Complex impedance method for AC circuits <https://physics.umd.edu/~jacobson/273c/impedance.pdf>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols

current = symbols.current
"""
Phasor :symbols:`current` through the component.
"""

voltage = symbols.voltage
"""
Phasor :symbols:`voltage` across the component.
"""

impedance = symbols.electrical_impedance
"""
:symbols:`electrical_impedance` of the component.
"""

law = Eq(current, voltage / impedance)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(voltage_=voltage, impedance_=impedance)
@validate_output(current)
def calculate_current(voltage_: Quantity, impedance_: Quantity) -> Quantity:
    result = solve(law, current)[0].subs({
        voltage: voltage_,
        impedance: impedance_,
    })
    return Quantity(result)


# UNIQUE_LAW_ID: 487
