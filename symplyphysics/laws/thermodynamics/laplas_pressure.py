from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)


surface_tension_of_the_liquid = Symbol("surface_tension_of_the_liquid", units.force / units.length)
radius_of_curvature = Symbol("radius_of_curvature", units.length)
laplas_pressure = Symbol("laplas_pressure", units.pressure)


law = Eq(laplas_pressure, 2 * surface_tension_of_the_liquid / radius_of_curvature)


def print_law() -> str:
    return print_expression(law)


@validate_input(surface_tension_of_the_liquid_=surface_tension_of_the_liquid,
                radius_of_curvature_=radius_of_curvature)
@validate_output(laplas_pressure)
def calculate_laplas_pressure(surface_tension_of_the_liquid_: Quantity,
                              radius_of_curvature_: Quantity) -> Quantity:
    solved = solve(law, laplas_pressure, dict=True)[0][laplas_pressure]
    result_expr = solved.subs({
        surface_tension_of_the_liquid: surface_tension_of_the_liquid_,
        radius_of_curvature: radius_of_curvature_
    })
    return Quantity(result_expr)
