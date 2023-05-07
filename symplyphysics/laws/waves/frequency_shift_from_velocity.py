from sympy import (Eq, solve, sqrt, simplify)
from sympy.physics.units import speed_of_light
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.relativistic.waves import longitudinal_frequency_shift_from_velocity as relativistic_doppler

# Description
## The Doppler effect or Doppler shift is the apparent change in frequency of a wave in relation to an observer moving relative to the wave source.

# Law: fo = fs * (v + vo)/(v + vs), where
## fo is observed frequency,
## fs is source wave frequency,
## v is wave velocity in this media,
## vo is observer velocity relative to the medium (positive when moving towards source, negative when moving away from source),
## vs is source velocity relative to the medium (positive when moving away from observer, negative when moving towards observer).

# Conditions:
## - Velocities vs and vo are much smaller than v.
## - Source and observer are moving directly towards or away from each other.

# Note:
## When direction of source and observer is not directed towards to/away from each other, one should introduce angle to the formula.

observed_frequency = Symbol("observed_frequency", units.frequency)
real_frequency = Symbol("real_frequency", units.frequency)
wave_velocity = Symbol("wave_velocity", units.velocity)
source_velocity = Symbol("source_velocity", units.velocity)
observer_velocity = Symbol("observer_velocity", units.velocity)

law = Eq(observed_frequency,
    real_frequency * (wave_velocity + observer_velocity) / (wave_velocity + source_velocity))

# Confirm that classical Doppler effect is a special case of relativistic Doppler effect

# This is a general formula for Doppler effect that has both classical and relativistic parts of equation.
# Note, that observer velocity sign is minus, meaning that observer velocity is positive, when moving away from source.
#TODO: reuse in relativistic Doppler effect
general_doppler_law = Eq(observed_frequency,
    real_frequency * (1 - observer_velocity / wave_velocity) / (1 + source_velocity / wave_velocity) * sqrt(
            (1 - (source_velocity / speed_of_light)**2) / (1 - (observer_velocity / speed_of_light)**2)))
# Relativistic part is zero for velocities much less than speed of light
classical_law = general_doppler_law.subs({(source_velocity / speed_of_light)**2: 0, (observer_velocity / speed_of_light)**2: 0})
# Change observer velocity sign to conform to our velocity notation
assert expr_equals(classical_law.subs(observer_velocity, -observer_velocity).rhs, law.rhs)

# As the signal propagation speed wave_velocity goes to speed_of_light, Doppler effect
# formula evolves to relativistic version
#TODO: move this proof to longitudinal_frequency_shift_from_velocity law
relativistic_law = simplify(general_doppler_law.subs(wave_velocity, speed_of_light))
# Relative velocity is a relativistic version of velocities addition
relative_velocity = (observer_velocity + source_velocity) / (1 + observer_velocity * source_velocity / speed_of_light**2)
applied_relativistic_law = relativistic_doppler.law.rhs.subs({
    relativistic_doppler.real_frequency: real_frequency,
    relativistic_doppler.source_velocity: relative_velocity})
# We verify that expressions inside square root are identical - that's enough to prove
# that our relativistic version of law is indeed a special case of general_doppler_law
assert expr_equals((applied_relativistic_law / real_frequency)**2, (relativistic_law.rhs / real_frequency)**2)


def print() -> str:
    return print_expression(law)


@validate_input_symbols(real_frequency_=real_frequency,
    wave_velocity_=wave_velocity,
    source_velocity_=source_velocity,
    observer_velocity_=observer_velocity)
@validate_output_symbol(observed_frequency)
def calculate_observed_frequency(real_frequency_: Quantity, wave_velocity_: Quantity,
    source_velocity_: Quantity, observer_velocity_: Quantity) -> Quantity:
    result_expr = solve(law, observed_frequency, dict=True)[0][observed_frequency]
    frequency_applied = result_expr.subs({
        real_frequency: real_frequency_,
        wave_velocity: wave_velocity_,
        source_velocity: source_velocity_,
        observer_velocity: observer_velocity_
    })
    return expr_to_quantity(frequency_applied)
