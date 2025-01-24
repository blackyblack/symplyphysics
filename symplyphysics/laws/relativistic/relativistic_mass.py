"""
Relativistic mass via rest mass and speed
=========================================

Relativistic mass is proportional to the rest mass of the object and approaches
infinity as its speed approaches the speed of light.

**Notation:**

#. :quantity_notation:`speed_of_light`.

**Conditions:**

#. Non-zero rest mass.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Mass_in_special_relativity#Relativistic_mass>`__.

..
    TODO rename file
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (Quantity, validate_input,
    validate_output, symbols, quantities)

relativistic_mass = symbols.relativistic_mass
"""
:symbols:`relativistic_mass` of the object.
"""

rest_mass = symbols.rest_mass
"""
:symbols:`rest_mass` of the object.
"""

speed = symbols.speed
"""
:symbols:`speed` of the object.
"""

law = Eq(relativistic_mass, rest_mass / sqrt(1 - speed**2 / quantities.speed_of_light**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(rest_mass_=rest_mass, velocity_=speed)
@validate_output(relativistic_mass)
def calculate_relativistic_mass(rest_mass_: Quantity, velocity_: Quantity) -> Quantity:
    result_expr = solve(law, relativistic_mass, dict=True)[0][relativistic_mass]
    mass_applied = result_expr.subs({rest_mass: rest_mass_, speed: velocity_})
    return Quantity(mass_applied)
