"""
Surface tension force via surface tension and length
====================================================

The surface tension force is directed tangentially to the surface of the fluid, perpendicular to
the section of the fluid contour on which it acts and is proportional to the length of this section.
Also see `figure <https://www.researchgate.net/publication/312093145/figure/fig7/AS:655095716921345@1533198398154/Illustration-of-surface-tension-as-a-force-per-unit-length-An-operator-extends-the-area.png>`_.
"""

from sympy import (Eq, solve)
from symplyphysics import (
    clone_symbol,
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

surface_tension_force = clone_symbol(symbols.dynamics.force)
"""
Surface tension :attr:`~symplyphysics.symbols.dynamics.force`.
"""

surface_tension = Symbol("surface_tension", units.force / units.length)
r"""
Surface tension of the fluid.

Symbol:
    :code:`gamma`

Latex:
    :math:`\gamma`
"""

length = Symbol("length", units.length)
"""
Length of the contour to which the force is applied.

Symbol:
    :code:`l`
"""

law = Eq(surface_tension_force, surface_tension * length)
r"""
:code:`F = gamma * l`

Latex:
    .. math::
        F = \gamma l
"""


@validate_input(surface_coefficient_=surface_tension, contour_length_=length)
@validate_output(surface_tension_force)
def calculate_force(surface_coefficient_: Quantity, contour_length_: Quantity) -> Quantity:
    result_expr = solve(law, surface_tension_force, dict=True)[0][surface_tension_force]
    result_expr = result_expr.subs({
        surface_tension: surface_coefficient_,
        length: contour_length_,
    })
    return Quantity(result_expr)
