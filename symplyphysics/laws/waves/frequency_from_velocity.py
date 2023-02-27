from symplyphysics import (
    symbols, Eq, pretty, solve, units, Quantity,
    validate_input, validate_output, expr_to_quantity
)

# Description
## If wave source moves in media, observed wave frequency differs from real. 
## If observer moves in media, observed wave frequency differs from real as well. This effect is also known as irrelative Dopler effect.

# Law: fo = fs * (v + vo)/(v + vs), where
## fo is observed frequency,
## fs is source wave frequency,
## v is wave velocity in this media,
## vo is observer velocity related to media,
## vs is source velocity related to media.

# Conditions:
## 1. velocities vs and vo are much smaller than v.
## 2. Positive velocity direction is observer to source.

observed_frequency = symbols('observed_frequency')
real_frequency = symbols('real_frequency')
wave_velocity = symbols('wave_velocity')
source_velocity = symbols('source_velocity')
observer_velocity = symbols('observer_velocity')

law = Eq(observed_frequency, real_frequency * (wave_velocity + observer_velocity) / (wave_velocity + source_velocity))

def print():
    return pretty(law, use_unicode=False)

@validate_input(real_freq_=units.frequency, wave_vel_=units.velocity, source_vel_=units.velocity, observer_vel_=units.velocity)
@validate_output(units.frequency)
def calculate_frequency(real_freq_: Quantity, wave_vel_: Quantity, source_vel_: Quantity, observer_vel_: Quantity) -> Quantity:        
    result_expr = solve(law, observed_frequency, dict=True)[0][observed_frequency]
    frequency_applied = result_expr.subs({real_frequency: real_freq_, 
                                          wave_velocity: wave_vel_,
                                          source_velocity: source_vel_,
                                          observer_velocity: observer_vel_})
    return expr_to_quantity(frequency_applied, 'observed_frequency')
