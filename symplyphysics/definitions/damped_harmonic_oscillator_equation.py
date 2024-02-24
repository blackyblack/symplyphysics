from sympy import Derivative, dsolve, solve, Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
    angle_type,
    dimensionless,
)

# Description
## Assuming there is a damping force acting on an oscillating body that is linearly proportional
## to the body's velocity, we can write a differential equation for the body's position. We're
## assuming the body only moves in one direction.

# Definition: d**2(x(t))/dt**2 + 2*zeta*omega*d(x(t))/dt + (omega**2)*x(t) = 0
## x(t) - position of the oscillating body
## t - time
## omega - undamped angular frequency
## zeta - damping ratio, the value of which critically determines the behavior of the system

displacement = Function("displacement", units.length)
time = Symbol("time", units.time, positive=True)
undamped_angular_frequency = Symbol("undamped_angular_frequency", angle_type / units.time, positive=True)
damping_ratio = Symbol("damping_ratio", dimensionless, positive=True)

definition = (
    Derivative(displacement(time), time, 2)
    + 2 * damping_ratio * undamped_angular_frequency * Derivative(displacement(time), time)
    + undamped_angular_frequency**2 * displacement(time)
)

# This solution can be used in the case of an overdamped oscillator.
# It contains coefficients "C1" and "C2" that can be found from initial conditions.
general_solution = dsolve(definition, displacement(time)).rhs


def print_law() -> str:
    return print_expression(definition)


@validate_input(
    initial_position_=units.length,
    initial_velocity_=units.velocity,
    undamped_angular_frequency_=undamped_angular_frequency,
    damping_ratio_=damping_ratio,
    time_=time,
)
@validate_output(displacement)
def calculate_displacement(
    initial_position_: Quantity,
    initial_velocity_: Quantity,
    undamped_angular_frequency_: Quantity,
    damping_ratio_: float,
    time_: Quantity,
) -> Quantity:
    dsolved = general_solution.subs({
        undamped_angular_frequency: undamped_angular_frequency_,
        damping_ratio: damping_ratio_,
    })
    c12 = solve(
        [
            Eq(initial_position_, dsolved.subs(time, 0)),
            Eq(initial_velocity_, dsolved.diff(time).subs(time, 0)),
        ],
        ("C1", "C2"),
    )
    dsolved = dsolved.subs(c12)
    result = dsolved.subs(time, time_)
    return Quantity(result)
