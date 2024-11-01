"""
Radiant exitance is radiant flux emitted per unit area
======================================================

**Radiant exitance**, or **radiant emittance** is the radiant flux emitted by a surface per unit area.

**Notation:**

#. The subscript :math:`e` stands for *energetic* to avoid confusion with photometric quantities.
"""

from sympy import Eq, Derivative
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
)

radiant_exitance = symbols.radiant_exitance
"""
:symbols:`radiant_exitance` of the surface.
"""

area = symbols.area
"""
:symbols:`area` of the surface.
"""

radiant_flux = clone_as_function(symbols.radiant_flux, [area])
"""
:symbols:`radiant_flux` emitted from the surface as a function of :attr:`~area`.
"""

definition = Eq(radiant_exitance, Derivative(radiant_flux(area), area))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    radiant_flux_=radiant_flux,
    surface_area_=area,
)
@validate_output(radiant_exitance)
def calculate_radiant_exitance(
    radiant_flux_: Quantity,
    surface_area_: Quantity,
) -> Quantity:
    # Calculate radiant exitance in case of constant radiant flux

    flux_function = radiant_flux_ * area / surface_area_
    result = definition.rhs.subs(radiant_flux(area), flux_function).doit()
    return Quantity(result)
