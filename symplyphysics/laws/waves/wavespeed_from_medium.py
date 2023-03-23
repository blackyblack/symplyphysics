from symplyphysics import (
    symbols, Eq, pretty, solve, units, Quantity, validate_output, expr_to_quantity
)

from sympy.physics.units import speed_of_light as C

# Description
## Wavespeed differs in different medium. Electromagnetic wave propagation speed depends on refraction factor of medium.
## Commonly refraction factor also depends on wave frequency.

# Law: Vmedium = C / n, where
## Vmedium is speed of electromagnetic wave in medium,
## C is speed of light in vacuum (it is a fundamental constant),
## n is refraction factor of medium.

wave_speed_in_medium = symbols("wavespeed_in_medium")
refraction_factor = symbols("refraction_factor")

law = Eq(wave_speed_in_medium, C / refraction_factor)

def print():
    return pretty(law, use_unicode=False)

@validate_output(units.velocity)
def calculate_wavespeed(refraction_factor_: Quantity) -> Quantity:        
    result_expr = solve(law, wave_speed_in_medium, dict=True)[0][wave_speed_in_medium]    
    wavespeed_applied = result_expr.subs({refraction_factor: refraction_factor_})    
    return expr_to_quantity(wavespeed_applied, 'wavespeed_in_medium')
