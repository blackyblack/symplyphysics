from sympy import Eq, solve, sqrt, log, pi
from sympy.physics.units import magnetic_constant
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

## Description
## The microstrip line is a dielectric substrate on which a metal strip is applied.
## To create a short circuit based on a microstrip line, a metallized hole can be made in
## the line. Such an opening will have a certain inductance, which can be calculated.

## Law is: L = (mu0 / (2 * pi)) * (h * ln((h + sqrt(r^2 + h^2)) / r) + 1.5 * (r - sqrt(r^2 + h^2))), where
## L - hole inductance,
## h - the thickness of the microstrip line substrate,
## r - radius of the metallized hole,
## mu0 - magnetic constant.

inductance = Symbol("inductance", units.inductance)

thickness_of_substrate = Symbol("thickness_of_substrate", units.length)
radius_of_hole = Symbol("radius_of_hole", units.length)

expression_1 = sqrt(radius_of_hole**2 + thickness_of_substrate**2)
expression_2 = thickness_of_substrate * log((thickness_of_substrate + expression_1) / radius_of_hole)
expression_3 = 1.5 * (radius_of_hole - expression_1)

law = Eq(inductance, (magnetic_constant / (2 * pi)) * (expression_2 + expression_3))


def print_law() -> str:
    return print_expression(law)


@validate_input(thickness_of_substrate_=thickness_of_substrate,
    radius_of_hole_=radius_of_hole,)
@validate_output(inductance)
def calculate_inductance(thickness_of_substrate_: Quantity, radius_of_hole_: Quantity) -> Quantity:
    result_expr = solve(law, inductance, dict=True)[0][inductance]
    result_expr = result_expr.subs({
        thickness_of_substrate: thickness_of_substrate_,
        radius_of_hole: radius_of_hole_,
    })
    return Quantity(result_expr)
