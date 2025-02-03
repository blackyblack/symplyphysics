"""
Excess pressure under curved surface of bubble
==============================================

Under the curved surface of the liquid, in addition to the internal pressure, additional
pressure is created due to the curvature of the surface. For a bubble, this pressure is
created by two surfaces: the outer and the inner.

**Links:**

#. `Physics LibreTexts, formula 20.2.4 <https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Classical_Mechanics_(Tatum)/20%3A_Miscellaneous/20.02%3A_Surface_Tension/20.2.01%3A_Excess_Pressure_Inside_Drops_and_Bubbles>`__.

..
    TODO: fix file name
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols, clone_as_symbol
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import laplace_pressure_of_spherical_shapes as laplace_law

pressure_difference = clone_as_symbol(symbols.pressure, display_symbol="Delta(p)", display_latex="\\Delta p")
"""
:symbols:`pressure` difference between the surfaces of the bubble.
"""

surface_tension = symbols.surface_tension
"""
:symbols:`surface_tension` of the bubble.
"""

radius = symbols.radius
"""
:symbols:`radius` of the bubble.
"""

law = Eq(pressure_difference, 4 * surface_tension / radius)
"""
:laws:symbol::

:laws:latex::
"""

# This law might be derived via Laplace law.

# TODO prefix all variables used in proof with underscore

_laplace_law_applied = laplace_law.law.subs({
    laplace_law.surface_tension: surface_tension,
    laplace_law.radius_of_curvature: radius
})

# This is an excess pressure in a spherical drop.
_pressure_derived = solve(_laplace_law_applied, laplace_law.laplace_pressure,
    dict=True)[0][laplace_law.laplace_pressure]

# Check if derived pressure is same as declared.
# The bubble has two surfaces – an outer and an inner one, each of which creates additional pressure.
# Therefore, an additional multiplier of "2" appears.
assert expr_equals(2 * _pressure_derived, law.rhs)


@validate_input(surface_tension_of_the_liquid_=surface_tension,
    radius_of_bubble_=radius)
@validate_output(pressure_difference)
def calculate_excessive_pressure(surface_tension_of_the_liquid_: Quantity,
    radius_of_bubble_: Quantity) -> Quantity:
    solved = solve(law, pressure_difference, dict=True)[0][pressure_difference]
    result_expr = solved.subs({
        surface_tension: surface_tension_of_the_liquid_,
        radius: radius_of_bubble_
    })
    return Quantity(result_expr)
