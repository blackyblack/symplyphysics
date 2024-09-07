r"""
Electrostatic potential energy of two charges via distance
==========================================================

Electrostatic potential energy due to two point charges depends on the inverse 
distance to the distance between the charges. Note that this is the energy of
interaction belonging to the entire system.
"""

from sympy import (Eq, solve, pi)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

electrostatic_potential_energy = Symbol("electrostatic_potential_energy", units.energy)
"""
Electrostatic potential energy of system.

Symbol:
    :code:`U_E`

Latex:
    :math:`U_E`
"""

absolute_permittivity = Symbol("absolute_permittivity", units.capacitance / units.length)
r"""
Absolute permittivity of the medium.

Symbol:
    :code:`epsilon`

Latex:
    :math:`\varepsilon`
"""

distance = Symbol("distance", units.length)
"""
Distance between the point charges.

Symbol:
    :code:`r`
"""

first_charge = Symbol("first_charge", units.charge)
"""
Value of the first charge.

Symbol:
    :code:`q1`

Latex:
    :math:`q_1`
"""

second_charge = Symbol("second_charge", units.charge)
"""
Value of the second charge.

Symbol:
    :code:`q2`

Latex:
    :math:`q_2`
"""

law = Eq(electrostatic_potential_energy,
    (first_charge * second_charge) / (4 * pi * absolute_permittivity * distance))
r"""
:code:`U_E = q1 * q2 / (4 * pi * epsilon * r)`

Latex:
    .. math::
        U_E = \frac{q_1 q_2}{4 \pi \varepsilon r}
"""

@validate_input(absolute_permittivity_=absolute_permittivity,
    distance_=distance,
    charge_1_=first_charge,
    charge_2_=second_charge)
@validate_output(electrostatic_potential_energy)
def calculate_energy(absolute_permittivity_: Quantity, distance_: Quantity, charge_1_: Quantity,
    charge_2_: Quantity) -> Quantity:
    result_expr = solve(law, electrostatic_potential_energy, dict=True)[0][electrostatic_potential_energy]
    result_expr = result_expr.subs({
        absolute_permittivity: absolute_permittivity_,
        distance: distance_,
        first_charge: charge_1_,
        second_charge: charge_2_,
    })
    return Quantity(result_expr)
