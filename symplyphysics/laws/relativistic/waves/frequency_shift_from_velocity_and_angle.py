import numbers
from sympy import (Eq, cos, solve, sqrt)
from sympy.physics.units import speed_of_light
from symplyphysics import (angle_type, units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input, validate_output)

# Description
## See [doppler effect](./longitudinal_frequency_shift_from_velocity.py) description. When objects are not moving collinear, one
## should add angles to the formula.

# Law: fo = fs * sqrt(c**2 - v**2)/(c - v * cos(phr)), where
## fo is observed frequency,
## fs is source wave frequency,
## v is relative speed of source and observer (magnitude of velocity),
## phr is angle between signal vector directed from source to observer, and source velocity,
## c is speed of light.

# Conditions:
## - Angle is detected at the moment of emission with respect to the observer frame.
## - Motion is in 2-D space.

# Note:
## It is not trivial to substitute moving source and idle observer with moving observer and idle source. Law
## depends on the current frame (observer or source are at rest at this frame), angle is detected at the point
## of emission or at the point of reception (see relativistic angle aberration).

observed_frequency = Symbol("observed_frequency", units.frequency)
real_frequency = Symbol("real_frequency", units.frequency)
relative_speed = Symbol("relative_speed", units.velocity)
source_angle = Symbol("source_angle", angle_type)

law = Eq(
    observed_frequency,
    real_frequency * sqrt(speed_of_light**2 - relative_speed**2) /
    (speed_of_light - relative_speed * cos(source_angle)))


def print() -> str:
    return print_expression(law)


@validate_input(real_frequency_=real_frequency,
    relative_speed_=relative_speed,
    source_angle_=source_angle)
@validate_output(observed_frequency)
def calculate_observed_frequency(real_frequency_: Quantity, relative_speed_: Quantity,
    source_angle_: float | Quantity) -> Quantity:
    #HACK: sympy angles are always in radians
    source_angle_radians = source_angle_ if isinstance(source_angle_,
        numbers.Number) else source_angle_.scale_factor
    result_expr = solve(law, observed_frequency, dict=True)[0][observed_frequency]
    frequency_applied = result_expr.subs({
        real_frequency: real_frequency_,
        relative_speed: relative_speed_,
        source_angle: source_angle_radians
    })
    return expr_to_quantity(frequency_applied)
