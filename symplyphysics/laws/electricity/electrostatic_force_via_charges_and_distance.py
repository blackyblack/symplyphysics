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
from symplyphysics import (clone_symbol, symbols, units, Quantity, Symbol,
    validate_input, validate_output)

electrostatic_force = clone_symbol(symbols.dynamics.force)
"""
Electrostatic force between two charges.

Symbol:
    :code:`F`
"""

first_charge = Symbol("first_charge", units.charge)
r"""
First charge.

Symbol:
    :code:`q_1`

Latex:
    :math:`q_1`
"""

second_charge = Symbol("second_charge", units.charge)
r"""
Second charge.

Symbol:
    :code:`q_2`

Latex:
    :math:`q_2`
"""

distance = Symbol("distance", units.length)
"""
Distance between the charges.

Symbol:
    :code:`r`
"""

law = Eq(electrostatic_force,
    1 / (4 * pi * units.vacuum_permittivity) * first_charge * second_charge / (distance**2))
r"""
:code:`F = 1 / (4 * pi * epsilon_0) * q_1 * q_2 / r^2`

Latex:
    .. math::
        F = \frac{1}{4 \pi \varepsilon_0} \frac{q_1 q_2}{r^2}
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
