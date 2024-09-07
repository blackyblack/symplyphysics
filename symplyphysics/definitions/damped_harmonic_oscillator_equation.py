"""
Damped harmonic oscillator equation
===================================

Assuming there is a damping force acting on an oscillating body that is linearly proportional
to the body's velocity, we can write a differential equation for the body's position. We're
assuming the body only moves in one direction.
"""

from sympy import Derivative, Eq, dsolve
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    FunctionNew,
    validate_input,
    validate_output,
    angle_type,
    dimensionless,
)

displacement = FunctionNew("x(t)", units.length, display_latex="x")
"""
Displacement of the oscillating body as a function of time.
"""

time = SymbolNew("t", units.time, positive=True)
"""
Time.
"""

undamped_angular_frequency = SymbolNew("w",
    angle_type / units.time,
    display_latex="\\omega",
    positive=True)
"""
Undamped angular frequency of the oscillator.
"""

damping_ratio = SymbolNew("z", dimensionless, display_latex="\\zeta", positive=True)
"""
Damping ratio, which critically determines the behavior of the system.
"""

definition = Eq(
    Derivative(displacement(time), time, 2) +
        2 * damping_ratio * undamped_angular_frequency * Derivative(displacement(time), time) +
        undamped_angular_frequency**2 * displacement(time),
    0)
"""
:laws:symbol::

:laws:latex::
"""


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
    ics = {
        displacement(0): initial_position_,
        displacement(time).diff(time).subs(time, 0): initial_velocity_,
    }
    dsolved = dsolve(
        definition,
        displacement(time),
        ics=ics,
    ).rhs
    result = dsolved.subs(time, time_).subs({
        undamped_angular_frequency: undamped_angular_frequency_,
        damping_ratio: damping_ratio_,
    })
    return Quantity(result)
