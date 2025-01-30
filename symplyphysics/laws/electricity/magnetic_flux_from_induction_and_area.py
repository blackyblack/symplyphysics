"""
Magnetic flux from magnetic flux density and area
=================================================

Magnetic flux is the flux of the vector of magnetic flux density through a certain surface.

**Links:**

#. `Wikipedia, first formula <https://en.wikipedia.org/wiki/Magnetic_flux>`__.
"""

from sympy import (Eq, solve, cos)
from symplyphysics import (Quantity, validate_input, validate_output, symbols)

magnetic_flux = symbols.magnetic_flux
"""
:symbols:`magnetic_flux` through the area.
"""

magnetic_flux_density = symbols.magnetic_flux_density
"""
:symbols:`magnetic_flux_density`.
"""

area = symbols.area
"""
:symbols:`area`.
"""

angle = symbols.angle
"""
:symbols:`angle` between the vector of magnetic flux density and the area vector.
"""

law = Eq(magnetic_flux, magnetic_flux_density * area * cos(angle))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(induction_=magnetic_flux_density, area_=area, angle_=angle)
@validate_output(magnetic_flux)
def calculate_flux(induction_: Quantity, area_: Quantity, angle_: Quantity | float) -> Quantity:
    result_flux_expr = solve(law, magnetic_flux, dict=True)[0][magnetic_flux]
    result_expr = result_flux_expr.subs({magnetic_flux_density: induction_, area: area_, angle: angle_})
    return Quantity(result_expr)
