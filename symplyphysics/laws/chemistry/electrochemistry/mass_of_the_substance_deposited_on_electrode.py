"""
Mass of substance deposited on electrode
========================================

Faraday's first law of electrolysis: the mass of a substance deposited on an electrode
during electrolysis is directly proportional to the amount of electricity transferred to
this electrode. By the amount of electricity, we mean the total electric charge that has
passed through the surface of the electrode.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Faraday%27s_laws_of_electrolysis#First_law>`__.

..
    TODO: replace `I * t` with charge `q`?
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols

mass = symbols.mass
"""
:symbols:`mass` deposited on the electrode.
"""

equivalent = symbols.electrochemical_equivalent
"""
:symbols:`electrochemical_equivalent`.
"""

current = symbols.current
"""
:symbols:`current`.
"""

time = symbols.time
"""
:symbols:`time`.
"""

law = Eq(mass, equivalent * current * time)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(equivalent_=equivalent, current_=current, time_=time)
@validate_output(mass)
def calculate_mass(equivalent_: Quantity, current_: Quantity, time_: Quantity) -> Quantity:
    result_expr = solve(law, mass, dict=True)[0][mass]
    result_expr = result_expr.subs({equivalent: equivalent_, current: current_, time: time_})
    return Quantity(result_expr)
