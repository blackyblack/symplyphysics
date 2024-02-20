from sympy import Derivative, Eq, symbols, Function as SymFunction
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
displacement = symbols("displacement", cls=SymFunction)
time = Symbol("time", units.time)
undamped_angular_frequency = Symbol("undamped_angular_frequency", angle_type / units.time)
damping_ratio = Symbol("damping_ratio", dimensionless)

definition = (
    Derivative(displacement(time), time, time)
    + 2 * damping_ratio * undamped_angular_frequency * Derivative(displacement(time), time)
    + undamped_angular_frequency**2 * displacement(time)
)


def print_law() -> str:
    return print_expression(definition)


# calculate_displacement?