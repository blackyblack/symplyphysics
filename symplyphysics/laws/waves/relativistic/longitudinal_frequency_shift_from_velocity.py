"""
Longitudinal frequency shift from speed
=======================================

Doppler effect is also applicable to electromagnetic waves in vacuum. As there is no any medium
required for these waves to propagate, speed of source related to observer is used for the
calculation of the Doppler effect.

**Notes:**

#. This is a special case, when wave speed in medium (vacuum) is close to speed of light. But
   relativistic version of classical Doppler effect is more general case, that includes
   non-electromagnetic waves.

**Conditions:**

#. Source and observer are moving directly towards or away from each other.
#. Wave speed is close to speed of light. It means this law is only applicable to
   electromagnetic waves.
#. Motion is in 1D space.
"""

from sympy import Eq, pi, solve, sqrt, simplify
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.quantities import speed_of_light
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.waves.relativistic import (
    longitudinal_frequency_shift_from_absolute_velocities as general_doppler_law,
)
from symplyphysics.laws.waves.relativistic import (
    frequency_shift_from_velocity_and_angle as relativistic_doppler_with_angle,
)

observer_frequency = clone_as_symbol(
    symbols.angular_frequency,
    display_symbol="f_o",
    display_latex="f_\\text{o}",
)
"""
Observer :symbols:`temporal_frequency`.
"""

source_frequency = clone_as_symbol(
    symbols.angular_frequency,
    display_symbol="f_s",
    display_latex="f_\\text{s}",
)
"""
Source :symbols:`temporal_frequency`.
"""

relative_speed = symbols.speed
"""
Relative speed between source and observer.
"""

law = Eq(
    observer_frequency,
    source_frequency * sqrt(
    (speed_of_light - relative_speed) / (speed_of_light + relative_speed)))
"""
:laws:symbol::

:laws:latex::
"""

# As the signal propagation speed wave_velocity goes to speed_of_light, Doppler effect
# formula evolves to relativistic version

_general_relativistic_law = simplify(
    general_doppler_law.law.subs({
    general_doppler_law.wave_speed: speed_of_light,
    general_doppler_law.source_frequency: source_frequency
    }))
# Relative velocity is a relativistic version of velocities addition
_add_velocities = (general_doppler_law.observer_speed +
    general_doppler_law.source_speed) / (1 +
    general_doppler_law.observer_speed * general_doppler_law.source_speed / speed_of_light**2)
_applied_law = law.rhs.subs(relative_speed, _add_velocities)
# We verify that expressions inside square root are identical - that's enough to prove
# that our relativistic version of law is indeed a special case of general_doppler_law
assert expr_equals((_general_relativistic_law.rhs / source_frequency)**2,
    (_applied_law / source_frequency)**2)

# Confirm that Doppler effect for collinear movement is a subset of Doppler effect with angles

## Classical Doppler effect angles are calculated with respect to the signal vector, directed
## from source to observer. Hence source moving directly towards observer has 0 angle. Therefore
## its cosine is positive, when moving towards observer.
## This law has reverse notation - source velocity is positive when moving away from the observer.
## Therefore we should use opposite direction for source - set pi as source angle.
_observed_frequency_zero_angles = relativistic_doppler_with_angle.law.rhs.subs({
    relativistic_doppler_with_angle.source_angle: pi,
    relativistic_doppler_with_angle.relative_speed: relative_speed,
    relativistic_doppler_with_angle.source_frequency: source_frequency,
})
## Square roots fail to compare with each other. Raise both parts to power of 2 before checking for equality.
assert expr_equals(_observed_frequency_zero_angles**2, law.rhs**2)


@validate_input(real_frequency_=source_frequency, relative_velocity_=relative_speed)
@validate_output(observer_frequency)
def calculate_observed_frequency(real_frequency_: Quantity,
    relative_velocity_: Quantity) -> Quantity:
    result_expr = solve(law, observer_frequency, dict=True)[0][observer_frequency]
    frequency_applied = result_expr.subs({
        source_frequency: real_frequency_,
        relative_speed: relative_velocity_
    })
    return Quantity(frequency_applied)
