r"""
Radiant exitance is radiant flux emitted per unit area
======================================================

*Radiant exitance*, or *radiant emittance* is the radiant flux emitted by a surface per unit area.

**Notation:**

#. The subscript :math:`e` stands for *energetical* to avoid confusion with photometric quantities.

#. :math:`\frac{\partial}{\partial A}` denotes a partial derivative w.r.t. area of the surface.
"""

from sympy import Eq, Derivative
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
)

radiant_exitance = Symbol("radiant_exitance", units.power / units.area)
"""
Radiant exitance of the surface.

Symbol:
    :code:`M_e`

Latex:
    :math:`M_e`
"""

radiant_flux = Function("radiant_flux", units.power)
r"""
Radiant flux, or radiant power, emitted from the surface.

Symbol:
    :code:`Phi_e`

Latex:
    :math:`\Phi_e`
"""

surface_area = Symbol("surface_area", units.area)
"""
The area of the surface.

Symbol:
    :code:`A`
"""

definition = Eq(radiant_exitance, Derivative(radiant_flux(surface_area), surface_area))
r"""
:code:`M_e = d(Phi_e)/dA`

Latex:
    .. math::
        M_e = \frac{\partial \Phi_e}{\partial A}
"""


@validate_input(
    radiant_flux_=radiant_flux,
    surface_area_=surface_area,
)
@validate_output(radiant_exitance)
def calculate_radiant_exitance(
    radiant_flux_: Quantity,
    surface_area_: Quantity,
) -> Quantity:
    # Calculate radiant exitance in case of constant radiant flux

    flux_function = radiant_flux_ * surface_area / surface_area_
    result = definition.rhs.subs(radiant_flux(surface_area), flux_function).doit()
    return Quantity(result)
