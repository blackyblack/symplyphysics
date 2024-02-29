from sympy import Eq, exp, cos, solve, symbols
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
    angle_type,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.symbols.quantities import scale_factor
from symplyphysics.definitions import damped_harmonic_oscillator_equation as damped_eqn
from symplyphysics.laws.kinematic.damped_oscillations import (
    damping_ratio_from_decay_constant_and_undamped_frequency as damping_ratio_law,
    damped_angular_frequency as damped_frequency_law,
)

# Description
## In the presence of a damping force in the oscillating system, the system's behaviour
## depends on the value of the damping ratio. When it is less than 1, the system oscillates
## with a slightly different frequency than in the undamped case, and its amplitude decreasing
## to zero. This behaviour is also known as underdamping.

# Law: x(t) = A * exp(-lambda*t) * cos(omega_damped*t + phi)
## x(t) - position of damped oscillator
## A - coefficient to be found from initial condition
## lambda - exponential decay constant
## omega_damped - [damped angular frequency](./damped_angular_frequency.py)
## phi - phase lag

# Conditions
## - System is underdamped, i.e. its damping ratio is less than 1.

displacement = Function("displacement", units.length, real=True)
time = Symbol("time", units.time, nonnegative=True)
coefficient = Symbol("coefficient", units.length, positive=True)
exponential_decay_constant = Symbol("exponential_decay_constant", 1 / units.time, positive=True)
damped_angular_frequency = Symbol("damped_angular_frequency",
    angle_type / units.time,
    positive=True)
phase_lag = Symbol("phase_lag", angle_type, real=True)

law = Eq(
    displacement(time),
    coefficient * exp(-1 * exponential_decay_constant * time) *
    cos(damped_angular_frequency * time + phase_lag))

# Relating to [damped oscillations equation](../../../definitions/damped_harmonic_oscillator_equation.py)
## We will show it is the solution of the aforementioned equation.

_damping_ratio_sym, _undamped_frequency_sym = symbols("damping_ratio undamped_frequency",
    positive=True)

_decay_eqn = damping_ratio_law.law.subs({
    damping_ratio_law.damping_ratio: _damping_ratio_sym,
    damping_ratio_law.exponential_decay_constant: exponential_decay_constant,
    damping_ratio_law.undamped_angular_frequency: _undamped_frequency_sym,
})
_frequencies_eqn = damped_frequency_law.law.subs({
    damped_frequency_law.damped_angular_frequency: damped_angular_frequency,
    damped_frequency_law.undamped_angular_frequency: _undamped_frequency_sym,
    damped_frequency_law.damping_ratio: _damping_ratio_sym,
})

_solved = solve([_decay_eqn, _frequencies_eqn], (_damping_ratio_sym, _undamped_frequency_sym),
    dict=True)[0]
_damping_ratio_expr = _solved[_damping_ratio_sym]
_undamped_frequency_expr = _solved[_undamped_frequency_sym]

_diff_eqn = damped_eqn.definition.subs(damped_eqn.time, time)

_diff_eqn_subs = _diff_eqn.subs({
    damped_eqn.displacement(time): law.rhs,
    damped_eqn.undamped_angular_frequency: _undamped_frequency_expr,
    damped_eqn.damping_ratio: _damping_ratio_expr,
}).doit()

assert expr_equals(_diff_eqn_subs, 0)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    initial_position_=units.length,
    exponential_decay_constant_=exponential_decay_constant,
    damped_angular_frequency_=damped_angular_frequency,
    phase_lag_=phase_lag,
    time_=time,
)
@validate_output(displacement)
def calculate_displacement(
    initial_position_: Quantity,
    exponential_decay_constant_: Quantity,
    damped_angular_frequency_: Quantity,
    phase_lag_: Quantity | float,
    time_: Quantity,
) -> Quantity:
    coefficient_expr = solve(
        law.subs({
        displacement(time): initial_position_,
        time: 0
        }),
        coefficient,
    )[0]
    displacement_expr = law.rhs.subs(coefficient, coefficient_expr)
    result = displacement_expr.subs({
        exponential_decay_constant: exponential_decay_constant_,
        damped_angular_frequency: damped_angular_frequency_,
        phase_lag: scale_factor(phase_lag_),
        time: time_,
    })
    return Quantity(result)
