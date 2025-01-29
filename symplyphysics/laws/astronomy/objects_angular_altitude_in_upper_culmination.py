"""
Object's angular altitude in upper culmination
==============================================

In observational astronomy, culmination is the passage of a celestial object across the
observer's local meridian. An object's angular altitude in degrees at its upper
culmination is equal to 90 minus the observer's latitude plus the object's declination.

**Conditions:**

#. The object is at its upper culmination.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Culmination#>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

altitude = symbols.altitude
"""
:symbols:`altitude` of the object.
"""

latitude = symbols.latitude
"""
:symbols:`latitude` of the object.
"""

declination = symbols.declination
"""
:symbols:`declination` of the object.
"""

ninety_degrees = Quantity(90 * units.deg, display_symbol="90_deg", display_latex="90^\\circ")
"""
A :math:`90^\\circ` angle.
"""

law = Eq(altitude, ninety_degrees - latitude + declination)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(latitude_=latitude, declination_=declination)
@validate_output(altitude)
def calculate_angular_altitude(latitude_: Quantity, declination_: Quantity) -> Quantity:
    result_expr = solve(law, altitude, dict=True)[0][altitude]
    result_expr = result_expr.subs({
        latitude: latitude_,
        declination: declination_,
    })
    return Quantity(result_expr)
