from sympy import (Eq, solve, S)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    dimensionless,
    convert_to,
    validate_input,
    validate_output,
)

# Description
## If wave transfers from one medium to another, it refracts. That's because of different propagation speeds in different mediums.
## This factor shows how much slower wave propagates in refracting medium related to outer medium.

# Definition: n = Vouter/Vrefraction
# Where:
## n is refractive index,
## Vouter is wave propagation speed in outer medium,
## Vrefracting is wave propagation speed in refracting medium.

# Conditions
## - Mediums are isotropic and transparent.
## - Wave is monochromic as propagation speed depends on frequency.

refractive_index = Symbol("refractive_index", dimensionless)
outer_speed = Symbol("outer_speed", units.velocity)
refracting_speed = Symbol("refracting_speed", units.velocity)

definition = Eq(refractive_index, outer_speed / refracting_speed)

definition_units_SI = S.One


def print_law() -> str:
    return print_expression(definition)


@validate_input(outer_speed_=outer_speed, refracting_speed_=refracting_speed)
@validate_output(refractive_index)
def calculate_refractive_index(outer_speed_: Quantity, refracting_speed_: Quantity) -> float:
    result_index_expr = solve(definition, refractive_index, dict=True)[0][refractive_index]
    result_expr = result_index_expr.subs({
        outer_speed: outer_speed_,
        refracting_speed: refracting_speed_
    })
    result = Quantity(result_expr)
    return float(convert_to(result, S.One).evalf())
