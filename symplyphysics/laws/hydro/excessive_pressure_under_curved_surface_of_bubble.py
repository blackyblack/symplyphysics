from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.laws.thermodynamics import laplas_pressure as laplas_law

# Description
## Under the curved surface of the liquid, in addition to the internal pressure,
## additional pressure is created due to the curvature of the surface.
## For a bubble, this pressure is created by two surfaces: the outer and the inner.

## Law is: p = 4 * sigma / R, where
## p - excessive pressure under the curved surface of bubble,
## sigma - surface tension of the liquid,
## R - radius of bubble.

excessive_pressure = Symbol("excessive_pressure", units.pressure)

surface_tension_of_the_liquid = Symbol("surface_tension_of_the_liquid", units.force / units.length)
radius_of_bubble = Symbol("radius_of_bubble", units.length)

law = Eq(excessive_pressure, 4 * surface_tension_of_the_liquid / radius_of_bubble)

# This law might be derived via Laplas law.
# The bubble has two surfaces – an outer and an inner one, each of which creates additional pressure.
# Therefore, an additional multiplier of "2" appears.

laplas_law_applied = laplas_law.law.subs({
    laplas_law.surface_tension_of_the_liquid: 2 * surface_tension_of_the_liquid,
    laplas_law.radius_of_curvature: radius_of_bubble
})
pressure_derived = solve(laplas_law_applied, laplas_law.laplas_pressure, dict=True)[0][laplas_law.laplas_pressure]

# Check if derived pressure is same as declared.
assert expr_equals(pressure_derived, law.rhs)


def print_law() -> str:
    return print_expression(law)


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
