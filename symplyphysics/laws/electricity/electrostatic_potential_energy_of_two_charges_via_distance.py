r"""
Electrostatic potential energy of two charges via distance
==========================================================

Electrostatic potential energy due to two point charges depends on the inverse 
distance to the distance between the charges. Note that this is the energy of
interaction belonging to the entire system.

**Notation:**

#. :math:`\varepsilon_0` (:code:`epsilon_0`) is vacuum permittivity.
"""

from sympy import (Eq, solve, pi)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    dimensionless,
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

relative_permittivity = Symbol("relative_permittivity", dimensionless)
r"""
Relative permittivity of the medium.

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
r"""
Value of the first charge.

Symbol:
    :code:`q_1`

Latex:
    :math:`q_1`
"""

second_charge = Symbol("second_charge", units.charge)
r"""
Value of the second charge.

Symbol:
    :code:`q_2`

Latex:
    :math:`q_2`
"""

law = Eq(electrostatic_potential_energy,
    1 / (4 * pi * units.vacuum_permittivity) * (first_charge * second_charge) / (relative_permittivity * distance))
r"""
:code:`V = 1 / (4 * pi * epsilon_0) * q_1 * q_2 / (epsilon * r)`

Latex:
    .. math::
        V = \frac{1}{4 \pi \varepsilon_0} \frac{q_1 q_2}{\varepsilon r}
"""

@validate_input(relative_permittivity_=relative_permittivity,
    distance_=distance,
    charge_1_=first_charge,
    charge_2_=second_charge)
@validate_output(electrostatic_potential_energy)
def calculate_energy(relative_permittivity_: float, distance_: Quantity, charge_1_: Quantity,
    charge_2_: Quantity) -> Quantity:
    result_expr = solve(law, electrostatic_potential_energy, dict=True)[0][electrostatic_potential_energy]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        distance: distance_,
        first_charge: charge_1_,
        second_charge: charge_2_,
    })
    return Quantity(result_expr)
