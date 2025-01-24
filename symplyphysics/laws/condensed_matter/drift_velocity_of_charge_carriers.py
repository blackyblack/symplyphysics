"""
Drift velocity of charge carriers
=================================

Drift velocity is the average velocity attained by charged particles, such as electrons,
in a material due to an electric field. In general, an electron in a conductor will
propagate randomly at the Fermi velocity, resulting in an average velocity of zero.
Applying an electric field adds to this random motion a small net flow in one direction;
this is the drift.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Drift_velocity#>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols

drift_velocity = symbols.drift_velocity
"""
:symbols:`drift_velocity` of charge carriers.
"""

mobility = symbols.mobility
"""
:symbols:`mobility` of charge carriers.
"""

electric_field_strength = symbols.electric_field_strength
"""
:symbols:`electric_field_strength`.
"""

law = Eq(drift_velocity, mobility * electric_field_strength)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(charge_carriers_mobility_=mobility,
    electric_intensity_=electric_field_strength)
@validate_output(drift_velocity)
def calculate_velocity(charge_carriers_mobility_: Quantity,
    electric_intensity_: Quantity) -> Quantity:
    result_expr = solve(law, drift_velocity, dict=True)[0][drift_velocity]
    result_expr = result_expr.subs({
        mobility: charge_carriers_mobility_,
        electric_field_strength: electric_intensity_,
    })
    return Quantity(result_expr)
