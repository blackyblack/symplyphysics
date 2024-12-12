"""
Electric dipole moment is charge times distance
===============================================

Electric dipole moment is a measure of the separation of positive and negative
electrical charges within a system, i.e. it is a measure of the system's overall
polarity. Also see :doc:`vector counterpart <laws.electricity.vector.electric_dipole_moment_is_charge_times_displacement>`.

**Conditions:**

#. The system can be viewed as composed of two point charges.
#. The point charges are equal by magnitude.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Electric_dipole_moment#Elementary_definition>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

electric_dipole_moment = symbols.electric_dipole_moment
"""
:symbols:`electric_dipole_moment` of the system.
"""

charge = symbols.charge
"""
Magnitude of one the two point :symbols:`charge`-s comprising the system.
"""

distance = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` between point charges.
"""

law = Eq(electric_dipole_moment, charge * distance)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(charge_=charge, distance_=distance)
@validate_output(electric_dipole_moment)
def calculate_electric_moment(charge_: Quantity, distance_: Quantity) -> Quantity:
    result_expr = solve(law, electric_dipole_moment, dict=True)[0][electric_dipole_moment]
    result_expr = result_expr.subs({
        charge: charge_,
        distance: distance_,
    })
    return Quantity(result_expr)
