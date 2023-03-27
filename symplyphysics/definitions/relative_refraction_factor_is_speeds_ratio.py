from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, convert_to
)

# Description
## If wave transfers from one medium to another, it refracts. That's because of different propagation speeds in different mediums.
## This factor shows how much slower wave propagates in refracting medium related to outer medium.
## Definition: n = Vouter/Vrefraction, where
## n is refraction factor,
## Vouter is wave propagation speed of the outer medium,
## Vrefracting is wave propagation speed of the refracting medium.

# Conditions
## Mediums are isotropic and transparent.
## Wave is monochromic as propagation speed depends on frequency.

refraction_factor, outer_speed, refracting_speed = symbols("refracting_factor, outer_speed, refracting_speed")
definition = Eq(refraction_factor, outer_speed / refracting_speed)

def print():
    return pretty(definition, use_unicode=False)

@validate_input(v1_=units.velocity, v2_=units.velocity)
def calculate_refraction_factor(v1_: Quantity, v2_: Quantity) -> float:
    result_factor_expr = solve(definition, refraction_factor, dict=True)[0][refraction_factor]
    v1_arg = convert_to(v1_, units.meter/units.second).subs({units.meter: 1, units.second: 1}).evalf(4)
    v2_arg = convert_to(v2_, units.meter/units.second).subs({units.meter: 1, units.second: 1}).evalf(4)
    result_expr = result_factor_expr.subs({outer_speed: v1_arg, refracting_speed: v2_arg})
    return result_expr