from sympy import Eq, exp, dsolve, solve
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
from symplyphysics.definitions import damped_harmonic_oscillator_equation as damped_eqn

# Description
## In the presence of a damping force in the oscillating system, the system's behaviour
## depends on the value of the damping ratio. When it is equal to 1, the system returns
## to equilibrium as quickly as possible without oscillating, but overshoot can occur if
## initial velocity is nonzero. This behaviour is also called critical damping.

# Law: x(t) = exp(-omega*t) * (x0 + (v0 + x0*omega)*t)
## x(t) - position of critically damped oscillator
## t - time
## omega - oscillator frequency in an undamped state
## x0 - oscillator's initial position
## v0 - oscillator's initial velocity

# Conditions
## - System is critically damped, i.e. its damping ratio is equal to 1.

displacement = Function("displacement", units.length, real=True)
time = Symbol("time", units.time, nonnegative=True)
undamped_angular_frequency = Symbol("undamped_angular_frequency",
    angle_type / units.time,
    positive=True)
initial_position = Symbol("initial_position", units.length, real=True)
initial_velocity = Symbol("initial_velocity", units.velocity, real=True)

law = Eq(
    displacement(time),
    exp(-1 * undamped_angular_frequency * time) * (initial_position +
    (initial_velocity + initial_position * undamped_angular_frequency) * time),
)

# Derive from damped oscillator equation

_eqn = damped_eqn.definition.subs(damped_eqn.time, time).subs({
    damped_eqn.displacement(time): displacement(time),
    damped_eqn.undamped_angular_frequency: undamped_angular_frequency,
    damped_eqn.damping_ratio: 1,
})

_dsolved = dsolve(_eqn, displacement(time)).rhs

_initial_position_eqn = Eq(initial_position, _dsolved.subs(time, 0))
_initial_velocity_eqn = Eq(initial_velocity, _dsolved.diff(time).subs(time, 0))

_c12 = solve(
    [_initial_position_eqn, _initial_velocity_eqn],
    ("C1", "C2"),
    dict=True,
)[0]

_dsolved = _dsolved.subs(_c12)

assert expr_equals(law.rhs, _dsolved)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    initial_position_=initial_position,
    initial_velocity_=initial_velocity,
    undamped_angular_frequency_=undamped_angular_frequency,
    time_=time,
)
@validate_output(displacement)
def calculate_displacement(
    initial_position_: Quantity,
    initial_velocity_: Quantity,
    undamped_angular_frequency_: Quantity,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        initial_position: initial_position_,
        initial_velocity: initial_velocity_,
        undamped_angular_frequency: undamped_angular_frequency_,
        time: time_,
    })
    return Quantity(result)
