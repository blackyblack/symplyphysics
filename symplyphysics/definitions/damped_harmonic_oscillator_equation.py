from sympy import (
    Derivative,
    symbols,
    Function as SymFunction,
    dsolve,
    sqrt,
    exp,
    cos,
)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    angle_type,
    dimensionless,
)
from symplyphysics.core.expr_comparisons import expr_equals

# Description
## Assuming there is a damping force acting on an oscillating body that is linearly proportional
## to the body's velocity, we can write a differential equation for the body's position. We're
## assuming the body only moves in one direction.

# Definition: d**2(x(t))/dt**2 + 2*zeta*omega*d(x(t))/dt + (omega**2)*x(t) = 0
## x(t) - position of the oscillating body
## t - time
## omega - undamped angular frequency
## zeta - damping ratio, the value of which critically determines the behavior of the system

# not specifying the dimensions because it can be any physical quantity
displacement = symbols("displacement", cls=SymFunction, real=True)
time = Symbol("time", units.time, positive=True)
undamped_angular_frequency = Symbol(
    "undamped_angular_frequency", angle_type / units.time, positive=True
)
damping_ratio = Symbol("damping_ratio", dimensionless, positive=True)

definition = (
    Derivative(displacement(time), time, time)
    + 2 * damping_ratio * undamped_angular_frequency * Derivative(displacement(time), time)
    + undamped_angular_frequency**2 * displacement(time)
)

# The general solution looks like this:
## x(t) = exp(-w*z*t) * (C1 * exp(w*sqrt(z**2 - 1)*t) + C2 * exp(-w*sqrt(z**2 - 1)*t))

# z > 1, overdamped system
## the system exponentially decays to steady state without oscillations
## larger values of the damping ratio return to equilibrium more slowly

overdamped_displacement = dsolve(definition, displacement(time)).rhs

# z = 1, critically damped system
## the system returns to equilibrium as quickly as possiblewithout oscillating
## but overshoot can occur if initial velocity is nonzero

critical_initial_position = symbols("critical_initial_position", real=True)
critical_initial_velocity = symbols("critical_initial_velocity", real=True)

# In this case the solution is in the form
## x_critical(t) = exp(-w*t) * (x0 + (v0 + x0*w)*t) where
## - x0 is the initial position of the oscillator
## - v0 is the initial velocity of the oscillator

critical_c1 = critical_initial_position
critical_c2 = critical_initial_velocity + critical_initial_position * undamped_angular_frequency

critical_displacement = (
    exp(-1 * undamped_angular_frequency * time)
    * (critical_c1 + critical_c2 * time)
)

_critical_subs = (
    definition
    .subs(damping_ratio, 1)
    .subs(displacement(time), critical_displacement)
    .doit()
    .simplify()
)
assert expr_equals(_critical_subs, 0)

# z < 1, underdamped
## the system oscillates with a slightly different frequency than the undamped case
## with the amplitude decreasing to zero

underdamped_magnitude = symbols("underdamped_magnitude", positive=True)
underdamped_angular_frequency = undamped_angular_frequency * sqrt(1 - damping_ratio**2)
underdamped_exponential_decay = undamped_angular_frequency * damping_ratio
underdamped_phase_lag = symbols("underdamped_phase_lag", real=True)

# The solition in the underdamped case take the form of
## x_underdamped(t) = x0 * exp(-lambda*t) * cos(w1*t + phi), where
## - x0 is the amplitude of underdamped oscillations,
## - lambda is the exponential decay,
## - phi is the phase lag

underdamped_displacement = (
    underdamped_magnitude
    * exp(-1 * underdamped_exponential_decay * time)
    * cos(underdamped_angular_frequency * time + underdamped_phase_lag)
)

_underdamped_subs = (
    definition
    .subs(displacement(time), underdamped_displacement)
    .doit()
    .simplify()
)
assert expr_equals(_underdamped_subs, 0)


def print_law() -> str:
    return print_expression(definition)


@validate_input(
    initial_position_=units.length,
    initial_velocity_=units.velocity,
    undamped_angular_frequency_=undamped_angular_frequency,
    time_=time,
)
@validate_output(units.length)
def calculate_critical_displacement(
    initial_position_: Quantity,
    initial_velocity_: Quantity,
    undamped_angular_frequency_: Quantity,
    damping_ratio_: float,
    time_: Quantity,
) -> Quantity:
    result = critical_displacement.subs({
        critical_initial_position: initial_position_,
        critical_initial_velocity: initial_velocity_,
        undamped_angular_frequency: undamped_angular_frequency_,
        damping_ratio: damping_ratio_,
    }).subs(time, time_)
    return Quantity(result)
