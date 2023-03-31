from symplyphysics import (
    symbols, Eq, pretty, Quantity, units, solve,
    validate_input, validate_output, expr_to_quantity, SI
)

# Description
## Wavelength is the spatial period of a periodic wave—the distance over which the wave's shape repeats.
## It is the distance between consecutive corresponding points of the same phase on the wave.
## Definition: λ = v * T, where
## λ is wavelength,
## v is wave spreading velocity (phase speed),
## T is oscillation period.

#TODO derive this from velocity and period definitions

wavelength = symbols('wavelength')
spreading_velocity, oscillation_period = symbols('spreading_velocity oscillation_period')

definition = Eq(wavelength, spreading_velocity * oscillation_period)

definition_dimension_SI = units.meter

def print():
    return pretty(definition, use_unicode=False)

def print_dimension():
    return pretty(definition_dimension_SI, use_unicode=False)

@validate_input(velocity_=units.velocity, period_=units.time)
@validate_output(units.length)
def calculate_wavelength(velocity_: Quantity, period_: Quantity) -> Quantity:    
    applied_definition = solve(definition, wavelength, dict=True)[0][wavelength]
    result_expr = applied_definition.subs({spreading_velocity: velocity_, oscillation_period: period_})    
    return expr_to_quantity(result_expr, 'calculated_wavelength')
