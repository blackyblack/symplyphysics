"""
Surface tension force via surface tension and length
====================================================

The surface tension force is directed tangentially to the surface of the fluid, perpendicular to
the section of the fluid contour on which it acts and is proportional to the length of this section.
Also see `figure <https://www.researchgate.net/publication/312093145/figure/fig7/AS:655095716921345@1533198398154/Illustration-of-surface-tension-as-a-force-per-unit-length-An-operator-extends-the-area.png>`_.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Surface_tension#Physics>`__.
"""

from sympy import Eq, solve
from symplyphysics import symbols, Quantity, validate_input, validate_output

surface_tension_force = symbols.force
"""
Surface tension :symbols:`force`.
"""

surface_tension = symbols.surface_tension
"""
:symbols:`surface_tension` of the fluid.
"""

length = symbols.length
"""
:symbols:`length` of the contour to which the force is applied.
"""

law = Eq(surface_tension_force, surface_tension * length)
"""
:laws:symbol::

:laws:latex::
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
