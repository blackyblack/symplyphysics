from symplyphysics import (
    symbols, Eq, pretty, solve, units, Quantity, validate_output, expr_to_quantity
)

# Description
## Wavespeed differs in different medium. Electromagnetic wave propagation speed depends on refraction factor of medium.

# Law: Vmedium = Vvacuum / n, where
## Vmedium is speed of electromagnetical wave in media,
## Vvacuum is speed of light in vacuum (it is a fundamentional constant),
## n is refraction factor of media.

# Conditions
## Commonly n depends on wave frequency.

wave_speed_in_media = symbols("wavespeed_in_media")
refraction_factor = symbols("refraction_factor")
c = symbols("speed_of_light")

c = 299792458 * units.meter / units.second

law = Eq(wave_speed_in_media, c / refraction_factor)

def print():
    return pretty(law, use_unicode=False)

@validate_output(units.velocity)
def calculate_wavespeed(refraction_factor_: Quantity) -> Quantity:        
    result_expr = solve(law, wave_speed_in_media, dict=True)[0][wave_speed_in_media]    
    wavespeed_applied = result_expr.subs({refraction_factor: refraction_factor_})    
    return expr_to_quantity(wavespeed_applied, 'wavespeed_in_media')
