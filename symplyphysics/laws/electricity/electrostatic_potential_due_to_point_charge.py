r"""
Electrostatic potential due to point charge
===========================================

Electrostatic potential of electric field due to a point charge is inversely proportional
to the distance to the point charge. Also see :doc:`laws.electricity.electrostatic_potential_is_work_to_bring_from_reference_point_over_charge`.

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

electrostatic_potential = Symbol("electrostatic_potential", units.voltage)
"""
Electrostatic potential at given point.

Symbol:
    :code:`V`
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
Distance to the point charge.

Symbol:
    :code:`r`
"""

charge = Symbol("charge", units.charge)
"""
Electric charge.

Symbol:
    :code:`q`
"""

law = Eq(electrostatic_potential,
    1 / (4 * pi * units.vacuum_permittivity) * charge / (relative_permittivity * distance))
r"""
:code:`V = 1 / (4 * pi * epsilon_0) * q / (epsilon * r)`

Latex:
    .. math::
        V = \frac{1}{4 \pi \varepsilon_0} \frac{q}{\varepsilon r}
"""


@validate_input(relative_permittivity_=relative_permittivity, distance_=distance, charge_=charge)
@validate_output(electrostatic_potential)
def calculate_electrostatic_potential(relative_permittivity_: float, distance_: Quantity,
    charge_: Quantity) -> Quantity:
    result_expr = solve(law, electrostatic_potential, dict=True)[0][electrostatic_potential]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        distance: distance_,
        charge: charge_,
    })
    return Quantity(result_expr)
