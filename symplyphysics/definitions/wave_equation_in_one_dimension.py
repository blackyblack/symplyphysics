from sympy import (
    Derivative,
    Eq,
    symbols,
    Function as SymFunction,
    cos,
)
from symplyphysics import (
    units,
    Symbol,
    print_expression,
)

# Description
## The wave equation is a second-order linear partial differential equation used to
## describe the propagation of waves, including standing wave fields such as mechanical
## or electromagnetic waves.

# Law: d**2(u(x, t))/dx**2 = (1/v**2) * d**2(u(x, t))/dt**2
## u(x, t) - factor representing a displacement from rest situation, 
##           which can be pressure, position, electric field, etc
## x - position
## t - time
## d**2/dx**2 - partial derivative w.r.t. position x
## d**2/dt**2 - partial derivative w.r.t. time t
## v - [phase velocity of wave](../laws/waves/phase_velocity_from_angular_frequency_and_wavenumber.py)

# Notes
## - This equation is called one-dimensional because the displacement function depends
##   only on one spatial dimension.

displacement = symbols("displacement", cls=SymFunction, real=True)
position = Symbol("position", units.length, real=True)
time = Symbol("time", units.time, positive=True)
phase_velocity = Symbol("phase_velocity", units.velocity, real=True)

definition = Eq(
    Derivative(displacement(position, time), position, 2),
    Derivative(displacement(position, time), time, 2) / phase_velocity**2,
)


def print_law() -> str:
    return print_expression(definition)
