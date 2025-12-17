"""
Electrostatic potential due to point charge
===========================================

Electrostatic potential of electric field due to a point charge is inversely proportional
to the distance to the point charge. Also see :doc:`laws.electricity.electrostatic_potential_is_work_to_bring_from_reference_point_over_charge`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Electric_potential#Electric_potential_due_to_a_point_charge>`__.
"""

from sympy import (Eq, solve, pi)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

electrostatic_potential = symbols.electric_potential
"""
Electrostatic potential at given point. See :symbols:`electric_potential`.
"""

absolute_permittivity = symbols.absolute_permittivity
"""
:symbols:`absolute_permittivity` of the medium.
"""

distance = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` to the point charge.
"""

charge = symbols.charge
"""
Electric :symbols:`charge`.
"""

law = Eq(electrostatic_potential, charge / (4 * pi * absolute_permittivity * distance))
"""
:laws:symbol::

:laws:latex::
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
