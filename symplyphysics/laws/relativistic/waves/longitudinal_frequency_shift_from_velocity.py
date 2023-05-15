from sympy import (Eq, solve, sqrt, simplify)
from sympy.physics.units import speed_of_light
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.relativistic.waves import longitudinal_frequency_shift_from_absolute_velocities as general_doppler_law

# Description
## Doppler effect is also applicable to electromagnetic waves in vacuum. As there is no any medium required for these waves to propagate,
## speed of source related to observer is used for the calculation of the Doppler effect.

# Law: fo = fs * sqrt((c - v)/(c + v)), where
## fo is observed frequency,
## fs is source wave frequency,
## v is relative velocity of source and observer (positive when moving away from each other, negative when moving towards each other),
## c is speed of light.

# Conditions:
## - Source and observer are moving directly towards or away from each other.
## - Wave speed is close to speed of light. It means this law is only applicable to electromagnetic waves.

# Note:
## When direction of source and observer is not directed towards to/away from each other, one should introduce angle to the formula to
## obtain transverse Doppler effect (TDE) formula.

# Note:
## This is a special case, when wave speed in medium (vacuum) is close to speed of light. But relativistic version of classical
## Doppler effect is more general case, that includes non-electromagnetic waves.

observed_frequency = Symbol("observed_frequency", units.frequency)
real_frequency = Symbol("real_frequency", units.frequency)
relative_velocity = Symbol("relative_velocity", units.velocity)

law = Eq(
    observed_frequency,
    real_frequency * sqrt((speed_of_light - relative_velocity) / (speed_of_light + relative_velocity)))

# As the signal propagation speed wave_velocity goes to speed_of_light, Doppler effect
# formula evolves to relativistic version

general_relativistic_law = simplify(general_doppler_law.law.subs({
    general_doppler_law.wave_velocity: speed_of_light,
    general_doppler_law.real_frequency: real_frequency}))
# Relative velocity is a relativistic version of velocities addition
add_velocities = (general_doppler_law.observer_velocity + general_doppler_law.source_velocity) / (1 + general_doppler_law.observer_velocity * general_doppler_law.source_velocity / speed_of_light**2)
applied_law = law.rhs.subs(relative_velocity, add_velocities)
# We verify that expressions inside square root are identical - that's enough to prove
# that our relativistic version of law is indeed a special case of general_doppler_law
assert expr_equals((general_relativistic_law.rhs / real_frequency)**2, (applied_law / real_frequency)**2)


def print() -> str:
    return print_expression(law)


@validate_input_symbols(real_frequency_=real_frequency, relative_velocity_=relative_velocity)
@validate_output_symbol(observed_frequency)
def calculate_observed_frequency(real_frequency_: Quantity, relative_velocity_: Quantity) -> Quantity:
    result_expr = solve(law, observed_frequency, dict=True)[0][observed_frequency]
    frequency_applied = result_expr.subs({
        real_frequency: real_frequency_,
        relative_velocity: relative_velocity_
    })
    return expr_to_quantity(frequency_applied)
