from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.laws.thermodynamics import laplace_pressure_of_spherical_shapes as laplace_law

# Description
## Under the curved surface of the liquid, in addition to the internal pressure,
## additional pressure is created due to the curvature of the surface.
## For a bubble, this pressure is created by two surfaces: the outer and the inner.

## Law is: p = 4 * sigma / R, where
## p - excessive pressure under the curved surface of bubble,
## sigma - surface tension of the liquid,
## R - radius of bubble.

# Links: Physics LibreTexts, formula 20.2.4 <https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Classical_Mechanics_(Tatum)/20%3A_Miscellaneous/20.02%3A_Surface_Tension/20.2.01%3A_Excess_Pressure_Inside_Drops_and_Bubbles>

excessive_pressure = Symbol("excessive_pressure", units.pressure)

surface_tension_of_the_liquid = Symbol("surface_tension_of_the_liquid", units.force / units.length)
radius_of_bubble = Symbol("radius_of_bubble", units.length)

law = Eq(excessive_pressure, 4 * surface_tension_of_the_liquid / radius_of_bubble)

# This law might be derived via Laplace law.

# TODO prefix all variables used in proof with underscore

laplace_law_applied = laplace_law.law.subs({
    laplace_law.surface_tension: surface_tension_of_the_liquid,
    laplace_law.radius_of_curvature: radius_of_bubble
})

# This is an excess pressure in a spherical drop.
pressure_derived = solve(laplace_law_applied, laplace_law.laplace_pressure,
    dict=True)[0][laplace_law.laplace_pressure]

# Check if derived pressure is same as declared.
# The bubble has two surfaces – an outer and an inner one, each of which creates additional pressure.
# Therefore, an additional multiplier of "2" appears.
assert expr_equals(2 * pressure_derived, law.rhs)


@validate_input(surface_tension_of_the_liquid_=surface_tension_of_the_liquid,
    radius_of_bubble_=radius_of_bubble)
@validate_output(excessive_pressure)
def calculate_excessive_pressure(surface_tension_of_the_liquid_: Quantity,
    radius_of_bubble_: Quantity) -> Quantity:
    solved = solve(law, excessive_pressure, dict=True)[0][excessive_pressure]
    result_expr = solved.subs({
        surface_tension_of_the_liquid: surface_tension_of_the_liquid_,
        radius_of_bubble: radius_of_bubble_
    })
    return Quantity(result_expr)
