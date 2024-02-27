from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## Under the curved surface of the liquid, in addition to the internal pressure, additional pressure is created due to the curvature of the surface.
## The excess pressure under the curved surface of the liquid is determined by the Laplace formula: P_l = 2 * sigma / R
## Where:
## P_l - Laplace pressure
## sigma - surface tension of the liquid
## R - radius of surface curvature

surface_tension_of_the_liquid = Symbol("surface_tension_of_the_liquid", units.force / units.length)
radius_of_curvature = Symbol("radius_of_curvature", units.length)
laplace_pressure = Symbol("laplace_pressure", units.pressure)

law = Eq(laplace_pressure, 2 * surface_tension_of_the_liquid / radius_of_curvature)


def print_law() -> str:
    return print_expression(law)


@validate_input(surface_tension_of_the_liquid_=surface_tension_of_the_liquid,
    radius_of_curvature_=radius_of_curvature)
@validate_output(laplace_pressure)
def calculate_laplace_pressure(surface_tension_of_the_liquid_: Quantity,
    radius_of_curvature_: Quantity) -> Quantity:
    solved = solve(law, laplace_pressure, dict=True)[0][laplace_pressure]
    result_expr = solved.subs({
        surface_tension_of_the_liquid: surface_tension_of_the_liquid_,
        radius_of_curvature: radius_of_curvature_
    })
    return Quantity(result_expr)
