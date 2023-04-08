from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, convert_to, S, expr_to_quantity
)

# Description
## If wave transfers from one medium to another, it refracts. That's because of different propagation speeds in different mediums.
## This factor shows how much slower wave propagates in refracting medium related to outer medium.
## Definition: n = Vouter/Vrefraction, where
## n is refractive index,
## Vouter is wave propagation speed in outer medium,
## Vrefracting is wave propagation speed in refracting medium.

# Conditions
## - Mediums are isotropic and transparent.
## - Wave is monochromic as propagation speed depends on frequency.

refractive_index, outer_speed, refracting_speed = symbols("refractive_index, outer_speed, refracting_speed")
definition = Eq(refractive_index, outer_speed / refracting_speed)

def print():
    return pretty(definition, use_unicode=False)

@validate_input(outer_speed_=units.velocity, refracting_speed_=units.velocity)
def calculate_refractive_index(outer_speed_: Quantity, refracting_speed_: Quantity) -> float:
    result_index_expr = solve(definition, refractive_index, dict=True)[0][refractive_index]    
    result_expr = result_index_expr.subs({outer_speed: outer_speed_, refracting_speed: refracting_speed_})
    return convert_to(expr_to_quantity(result_expr, 'refractive_index'), S.One)