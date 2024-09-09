"""
Current density via number density and drift velocity
=====================================================

Current density is the amount of charge per unit time that flows through a unit area of a chosen
cross section. The current density vector is defined as a vector whose magnitude is the electric
current per cross-sectional area at a given point in space, its direction being that of the motion
of the positive charges at this point.

**Conditions:**

#. Charge carriers carry the the same amount of charge.
"""

from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

current_density = Symbol("current_density", units.current / units.area)
"""
Current density of charge carriers.

Symbol:
    :code:`j`
"""

number_density = Symbol("number_density", 1 / units.volume)
"""
Number density, or number of charge carriers per unit volume.

Symbol:
    :code:`n`
"""

drift_velocity = Symbol("drift_velocity", units.velocity)
r"""
Drift velocity of charge carriers.

Symbol:
    :code:`v_d`

Latex:
    :math:`v_\text{d}`
"""

charge = Symbol("charge", units.charge)
"""
Charge of the charge carriers.

Symbol:
    :code:`q`
"""

law = Eq(current_density, charge * number_density * drift_velocity)
r"""
:code:`j = q * n * v_d`

Latex:
    .. math::
        j = q n v_\text{d}
"""


@validate_input(charge_carriers_concentration_=number_density,
    drift_velocity_=drift_velocity,
    charge_=charge)
@validate_output(current_density)
def calculate_current(charge_carriers_concentration_: Quantity, drift_velocity_: Quantity,
    charge_: Quantity) -> Quantity:
    result_expr = solve(law, current_density, dict=True)[0][current_density]
    result_expr = result_expr.subs({
        number_density: charge_carriers_concentration_,
        drift_velocity: drift_velocity_,
        charge: charge_
    })
    return Quantity(result_expr)
