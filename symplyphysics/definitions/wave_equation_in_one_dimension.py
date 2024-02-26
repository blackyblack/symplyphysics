from sympy import (
    Derivative,
    Eq,
    symbols,
    Function as SymFunction,
)
from symplyphysics import (
    units,
    Symbol,
    print_expression,
)

# Description
## The wave equation is a second-order linear partial differential equation used to
## describe the description of waves, including standing wave fields such as mechanical
## or electromagnetic waves.

# Law: d**2(u(x, t))/dx**2 = (1/v**2) * d**2(u(x, t))/dt**2
## u(x, t) - factor representing a displacement from rest situation
## x - position
## t - time
## d**2/dx**2 - partial derivative w.r.t. position x
## d**2/dt**2 - partial derivative w.r.t. time t
## v - wave speed (a fixed non-negative real value)

displacement = symbols("displacement", cls=SymFunction, real=True)
position = Symbol("position", units.length, real=True)
time = Symbol("time", units.time, positive=True)
wave_speed = Symbol("wave_speed", units.velocity, nonnegative=True)

definition = Eq(
    Derivative(displacement(position, time), position, 2),
    Derivative(displacement(position, time), time, 2) / wave_speed**2,
)


def print_law() -> str:
    return print_expression(definition)
