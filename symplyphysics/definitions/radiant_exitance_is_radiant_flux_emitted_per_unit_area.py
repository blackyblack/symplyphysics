"""
Radiant exitance is radiant flux emitted per unit area
======================================================

*Radiant exitance*, or *radiant emittance* is the radiant flux emitted by a surface per unit area.

**Notation:**

#. The subscript :math:`e` stands for *energetical* to avoid confusion with photometric quantities.
"""

from sympy import Eq, Derivative
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    FunctionNew,
    validate_input,
    validate_output,
)

radiant_exitance = SymbolNew("M_e", units.power / units.area)
"""
Radiant exitance of the surface.
"""

radiant_flux = FunctionNew("Phi_e(A)", units.power, display_latex="\\Phi_e")
"""
Radiant flux, or radiant power, emitted from the surface.
"""

area = SymbolNew("A", units.area)
"""
The area of the surface.
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
