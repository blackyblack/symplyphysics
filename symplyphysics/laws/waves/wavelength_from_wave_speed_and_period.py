from symplyphysics import (
    symbols, Eq, pretty, Quantity, units, solve,
    validate_input, validate_output, expr_to_quantity, SI
)

# Description
## Wavelength is the spatial period of a periodic wave — the distance over which the wave's shape repeats.
## It is the distance between consecutive corresponding points of the same phase on the wave.
## Law: λ = v * T, where
## λ is wavelength,
## v is wave propagation speed (phase speed),
## T is oscillation period.

#TODO derive this from velocity and period definitions

wavelength = symbols('wavelength')
propagation_speed, oscillation_period = symbols('propagation_speed oscillation_period')

law = Eq(wavelength, propagation_speed * oscillation_period)

def print():
    return pretty(law, use_unicode=False)

@validate_input(velocity_=units.velocity, period_=units.time)
@validate_output(units.length)
def calculate_wavelength(velocity_: Quantity, period_: Quantity) -> Quantity:    
    applied_definition = solve(law, wavelength, dict=True)[0][wavelength]
    result_expr = applied_definition.subs({propagation_speed: velocity_, oscillation_period: period_})    
    return expr_to_quantity(result_expr, 'calculated_wavelength')
