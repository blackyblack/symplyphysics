r"""
Electrostatic potential due to point charge
===========================================

Electrostatic potential of electric field due to a point charge is inversely proportional
to the distance to the point charge. Also see :doc:`laws.electricity.electrostatic_potential_is_work_to_bring_from_reference_point_over_charge`.
"""

from sympy import (Eq, solve, pi)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

electrostatic_potential = Symbol("electrostatic_potential", units.voltage)
"""
Electrostatic potential at given point.

Symbol:
    :code:`V`
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
    charge / (4 * pi * absolute_permittivity * distance))
r"""
:code:`V = q / (4 * pi * epsilon * r)`

Latex:
    .. math::
        V = \frac{q}{4 \pi \varepsilon r}
"""


@validate_input(absolute_permittivity_=absolute_permittivity, distance_=distance, charge_=charge)
@validate_output(electrostatic_potential)
def calculate_electrostatic_potential(absolute_permittivity_: Quantity, distance_: Quantity,
    charge_: Quantity) -> Quantity:
    result_expr = solve(law, electrostatic_potential, dict=True)[0][electrostatic_potential]
    result_expr = result_expr.subs({
        absolute_permittivity: absolute_permittivity_,
        distance: distance_,
        charge: charge_,
    })
    return Quantity(result_expr)
