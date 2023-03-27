from sympy import pi
from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## The main characteristics of any wave are it's frequency, wavelength, magnitude and propagation speed.
## Magnitude is power characteristic of wave.
## Other characteristics are dependent from each over.

## Definition: λ = V / (n * f)
## Where: λ is wavelength,
## V is propagation speed,
## n is refraction factor of the media,
## f is frequency.

wavelength, propagation_speed, refraction_factor, wave_frequency = symbols('wavelength, spreading_velocity, refraction_factor, wave_frequency')

definition = Eq(wavelength, propagation_speed / (refraction_factor * wave_frequency))

definition_dimension_SI = units.meter

def print():
    return pretty(definition, use_unicode=False)

def print_dimension():
    return pretty(definition_dimension_SI, use_unicode=False)

@validate_input(propagation_speed_=units.velocity, frequency_=units.frequency)
@validate_output(units.length)
def calculate_wavelength(propagation_speed_: Quantity, refraction_factor_: Quantity, frequency_: Quantity) -> Quantity:
    solved = solve(definition, wavelength, dict=True)[0][wavelength]
    result_expr = solved.subs({propagation_speed: propagation_speed_,
                               refraction_factor: refraction_factor_,
                               wave_frequency: frequency_})
    return expr_to_quantity(result_expr, 'wavelength')
