r"""
Electrostatic force via charges and distance
============================================

The *Coulomb's law* states that the electrostatic force between two point charges
in a vacuum is proportional to their values and inversely proportional to the
square of the distance between them.

**Notation:**

#. :math:`\varepsilon_0` (:code:`epsilon_0`) is vacuum permittivity.
"""

from sympy import (Eq, solve, pi)
from symplyphysics import (clone_as_symbol, symbols, Quantity,
    validate_input, validate_output, quantities)

electrostatic_force = symbols.force
"""
Electrostatic :symbols:`force` between two charges.
"""

first_charge = clone_as_symbol(symbols.charge, display_symbol="q_1")
"""
First charge.
"""

second_charge = clone_as_symbol(symbols.charge, display_symbol="q_2")
"""
Second charge.
"""

distance = symbols.distance_to_origin
"""
Distance between the charges.
"""

law = Eq(electrostatic_force,
    1 / (4 * pi * quantities.vacuum_permittivity) * first_charge * second_charge / (distance**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(first_charge_=first_charge, second_charge_=second_charge, distance_=distance)
@validate_output(electrostatic_force)
def calculate_force(first_charge_: Quantity, second_charge_: Quantity,
    distance_: Quantity) -> Quantity:
    solved = solve(law, electrostatic_force, dict=True)[0][electrostatic_force]
    result_expr = solved.subs({
        first_charge: first_charge_,
        second_charge: second_charge_,
        distance: distance_,
    })
    return Quantity(result_expr)
