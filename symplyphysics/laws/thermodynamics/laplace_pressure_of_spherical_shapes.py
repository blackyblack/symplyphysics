"""
Laplace pressure of spherical shapes
====================================

Under the curved surface of a fluid, in addition to the internal pressure, additional pressure
is created due to the curvature of the surface. The excess pressure under the curved surface of
the fluid is called *Laplace pressure*.

**Links:**

#. `Wikipedia, second formula <https://en.wikipedia.org/wiki/Laplace_pressure>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

laplace_pressure = clone_as_symbol(symbols.pressure, display_symbol="P_L", display_latex="P_\\text{L}")
"""
Excess :symbols:`pressure` under the curved surface.
"""

surface_tension = symbols.surface_tension
"""
:symbols:`surface_tension` of the fluid.
"""

radius_of_curvature = symbols.radius_of_curvature
"""
:symbols:`radius_of_curvature` of the surface.
"""

law = Eq(laplace_pressure, 2 * surface_tension / radius_of_curvature)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(surface_tension_of_the_liquid_=surface_tension,
    radius_of_curvature_=radius_of_curvature)
@validate_output(laplace_pressure)
def calculate_laplace_pressure(surface_tension_of_the_liquid_: Quantity,
    radius_of_curvature_: Quantity) -> Quantity:
    solved = solve(law, laplace_pressure, dict=True)[0][laplace_pressure]
    result_expr = solved.subs({
        surface_tension: surface_tension_of_the_liquid_,
        radius_of_curvature: radius_of_curvature_
    })
    return Quantity(result_expr)
