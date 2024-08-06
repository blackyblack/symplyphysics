"""
Laplace pressure of spherical shapes
====================================

Under the curved surface of a fluid, in addition to the internal pressure, additional pressure
is created due to the curvature of the surface. The excess pressure under the curved surface of
the fluid is called *Laplace pressure*.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

laplace_pressure = Symbol("laplace_pressure", units.pressure)
r"""
Excess pressure under the curved surface.

Symbol:
    :code:`P_L`

Latex:
    :math:`P_\text{L}`
"""

surface_tension = Symbol("surface_tension", units.force / units.length)
r"""
Surface tension of the fluid.

Symbol:
    :code:`sigma`

Latex:
    :math:`\sigma`
"""

radius_of_curvature = Symbol("radius_of_curvature", units.length)
"""
Radius of curvature of the surface.

Symbol:
    :code:`R`
"""

law = Eq(laplace_pressure, 2 * surface_tension / radius_of_curvature)
r"""
:code:`P_L = 2 * sigma / R`

Latex:
    .. math::
        P_\text{L} = \frac{2 \sigma}{R}
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
