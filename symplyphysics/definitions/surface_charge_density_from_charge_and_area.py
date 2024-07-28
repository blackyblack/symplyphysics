"""
Surface charge density from charge and area
===========================================

*Surface charge density* is the amount of charge per area unit.
It is a measure of how much electric charge is distributed throughout a surface.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

surface_charge_density = Symbol("surface_charge_density", units.charge / units.area)
r"""
Surface charge density of the surface.

Symbol:
    :code:`sigma`

Latex:
    :math:`\sigma`
"""

charge = Symbol("charge", units.charge)
"""
Total charge of the surface.

Symbol:
    :code:`q`
"""

area = Symbol("area", units.area)
"""
Area of the surface.

Symbol:
    :code:`A`
"""

definition = Eq(surface_charge_density, charge / area)
r"""
:code:`sigma = q / A`

Latex:
    .. math::
        \sigma = \frac{q}{A}
"""


@validate_input(charge_=charge, area_=area)
@validate_output(surface_charge_density)
def calculate_surface_charge_density(charge_: Quantity, area_: Quantity) -> Quantity:
    result_expr = solve(definition, surface_charge_density, dict=True)[0][surface_charge_density]
    result_surface_charge_density = result_expr.subs({
        charge: charge_,
        area: area_,
    })
    return Quantity(result_surface_charge_density)
