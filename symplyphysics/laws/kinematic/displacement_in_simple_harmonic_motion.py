from sympy import Eq, cos, symbols, Function as SymFunction, dsolve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    angle_type,
    validate_input,
    print_expression,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.quantity_decorator import validate_output_same
from symplyphysics.core.symbols.quantities import scale_factor
from symplyphysics.definitions import harmonic_oscillator_is_second_derivative_equation as harmonic_eqn

# Description
## Any motion that repeats at regular intervals is called periodic, or harmonic, motion.
## Simple harmonic motion is a particular type of repeated motion that is a sinusoidal
## function of time.

# Law: q(t) = q_m * cos(w * t + phi)
## q(t) - displacement function
## t - time
## q_m - oscillation amplitude
## w - angular frequency of oscillations
## phi - phase lag

displacement = symbols("displacement", cls=SymFunction)
amplitude = symbols("amplitude", positive=True)
angular_frequency = Symbol("angular_frequency", angle_type / units.time, positive=True)
phase_lag = Symbol("phase_lag", angle_type, real=True)
time = Symbol("time", units.time, real=True)

law = Eq(displacement(time), amplitude * cos(angular_frequency * time + phase_lag))

# Derive law from [oscillator equation](../../definitions/harmonic_oscillator_is_second_derivative_equation.py)

_eqn = harmonic_eqn.definition.subs({
    harmonic_eqn.time: time,
    harmonic_eqn.angular_frequency: angular_frequency,
}).subs(harmonic_eqn.displacement_function(time), displacement(time))

_initial_position = law.rhs.subs(time, 0)
_initial_velocity = law.rhs.diff(time).subs(time, 0)

_dsolved = dsolve(_eqn,
    displacement(time),
    ics={
    displacement(0): _initial_position,
    displacement(time).diff(time).subs(time, 0): _initial_velocity,
    }).rhs

assert expr_equals(law.rhs, _dsolved)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    angular_frequency_=angular_frequency,
    phase_lag_=phase_lag,
    time_=time,
)
@validate_output_same("amplitude_")
def calculate_displacement(
    amplitude_: Quantity,
    angular_frequency_: Quantity,
    phase_lag_: Quantity | float,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        amplitude: amplitude_,
        angular_frequency: angular_frequency_,
        phase_lag: scale_factor(phase_lag_),
        time: time_,
    })
    return Quantity(result)
