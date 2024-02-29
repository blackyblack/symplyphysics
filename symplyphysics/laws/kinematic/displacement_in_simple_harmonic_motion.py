from sympy import Eq, cos, sin, symbols, Function as SymFunction, dsolve, solve
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
}).subs(
    harmonic_eqn.displacement_function(time), displacement(time)
)

_A, _B = symbols("A B", real=True)

_dsolved = dsolve(_eqn, displacement(time)).rhs.subs({
    "C1": _A,
    "C2": _B,
})

# Setting amplitude equal to sqrt(A**2 + B**2) we can have these identities
_eqn_A = Eq(_A/amplitude, -1 * sin(phase_lag))
_eqn_B = Eq(_B/amplitude, cos(phase_lag))

_dsolved_subs = solve(
    [Eq(displacement(time), _dsolved), _eqn_A, _eqn_B],
    (displacement(time), _A, _B),
    dict=True
)[0][displacement(time)].simplify()

assert expr_equals(law.rhs, _dsolved_subs)


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
