import numbers
from sympy import (Eq, cos, solve)
from symplyphysics import (units, angle_type, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol)

# Description
## See [doppler effect](./frequency_shift_from_velocity.py) description. When objects are not moving collinear, one
## should add angles to the formula.

# Law: fo = fs * (v - vo * cos(pho))/(v - vs * cos(phs)), where
## fo is observed frequency,
## fs is source wave frequency,
## v is wave velocity in this media,
## vo is observer speed relative to the medium (magnitude of the velocity vector),
## vs is source speed relative to the medium (magnitude of the velocity vector),
## pho is angle between vector pointing from source of the wave to the observer (signal vector), and source velocity vector,
## phs is angle between signal vector and observer velocity vector.

# Conditions:
## - Source and observer velocities are less or equal than wave velocity. Otherwise emitted waves are left behind the source or never
## reach the observer.


observed_frequency = Symbol("observed_frequency", units.frequency)
real_frequency = Symbol("real_frequency", units.frequency)
wave_velocity = Symbol("wave_velocity", units.velocity)
source_speed = Symbol("source_speed", units.velocity)
observer_speed = Symbol("observer_speed", units.velocity)
source_angle = Symbol("source_angle", angle_type)
observer_angle = Symbol("observer_angle", angle_type)

law = Eq(observed_frequency,
    real_frequency * (wave_velocity - observer_speed * cos(observer_angle)) / (wave_velocity - source_speed * cos(source_angle)))

#TODO: add proof for Doppler effect


def print() -> str:
    return print_expression(law)


@validate_input_symbols(real_frequency_=real_frequency,
    wave_velocity_=wave_velocity,
    source_speed_=source_speed,
    observer_speed_=observer_speed,
    source_angle_=source_angle,
    observer_angle_=observer_angle)
@validate_output_symbol(observed_frequency)
def calculate_observed_frequency(real_frequency_: Quantity, wave_velocity_: Quantity,
    source_speed_: Quantity, observer_speed_: Quantity,
    source_angle_: float | Quantity, observer_angle_: float | Quantity) -> Quantity:
    #HACK: sympy angles are always in radians
    source_angle_radians = source_angle_ if isinstance(source_angle_, numbers.Number) else source_angle_.scale_factor
    observer_angle_radians = observer_angle_ if isinstance(observer_angle_, numbers.Number) else observer_angle_.scale_factor
    result_expr = solve(law, observed_frequency, dict=True)[0][observed_frequency]
    frequency_applied = result_expr.subs({
        real_frequency: real_frequency_,
        wave_velocity: wave_velocity_,
        source_speed: source_speed_,
        observer_speed: observer_speed_,
        source_angle: source_angle_radians,
        observer_angle: observer_angle_radians
    })
    return expr_to_quantity(frequency_applied)
