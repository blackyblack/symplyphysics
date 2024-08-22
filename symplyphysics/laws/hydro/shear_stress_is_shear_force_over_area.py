"""
Shear stress is shear force over area
=====================================

Shear stress is the component of stress coplanar with the material cross section on which it acts,
that is, it lies in the same plane as the area it is applied to. It is opposed to normal stress,
which arises from the force vector component perpendicular to the material cross section.
"""

from sympy import Eq, solve
from symplyphysics import (
    clone_symbol,
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

shear_stress = Symbol("shear_stress", units.pressure)
r"""
Shear stress.

Symbol:
    :code:`tau`

Latex:
    :math:`\tau`
"""

shear_force = clone_symbol(symbols.dynamics.force, "shear_force")
"""
Shear force applied.

Symbol:
    :code:`F`
"""

area = Symbol("area", units.area)
"""
Area of the material face parallel to the shear force applied.

Symbol:
    :code:`A`
"""

law = Eq(shear_stress, shear_force / area)
r"""
:code:`tau = F / A`

Latex:
    .. math::
        \tau = \frac{F}{A}
"""


@validate_input(force_applied_=shear_force, area_=area)
@validate_output(shear_stress)
def calculate_shear_stress(force_applied_: Quantity, area_: Quantity) -> Quantity:
    solved = solve(law, shear_stress)[0]
    result = solved.subs({
        shear_force: force_applied_,
        area: area_,
    })
    return Quantity(result)
