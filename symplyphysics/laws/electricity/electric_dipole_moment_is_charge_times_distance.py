"""
Electric dipole moment is charge times distance
===============================================

Electric dipole moment is a measure of the separation of positive and negative
electrical charges within a system, i.e. it is a measure of the system's overall
polarity. Also see :doc:`vector counterpart <laws.electricity.vector.electric_dipole_moment_is_charge_times_displacement>`.

**Conditions:**

#. The system can be viewed as composed of two point charges.
#. The point charges are equal by magnitude.
"""

from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

electric_dipole_moment = Symbol("electric_dipole_moment", units.charge * units.length)
"""
Electric dipole moment of the system.

Symbol:
    :code:`p`
"""

charge = Symbol("charge", units.charge)
"""
Magnitude of one the two point charges comprising the system.

Symbol:
    :code:`q`
"""

distance = Symbol("distance", units.length)
"""
Distance between point charges.

Symbol:
    :code:`d`
"""

law = Eq(electric_dipole_moment, charge * distance)
r"""
:code:`p = q * d`

Latex:
    .. math::
        p = q d
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
