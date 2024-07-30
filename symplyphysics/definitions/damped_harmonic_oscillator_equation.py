"""
Damped harmonic oscillator equation
===================================

Assuming there is a damping force acting on an oscillating body that is linearly proportional
to the body's velocity, we can write a differential equation for the body's position. We're
assuming the body only moves in one direction.
"""

from sympy import Derivative, dsolve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
    angle_type,
    dimensionless,
)

displacement = Function("displacement", units.length)
"""
Displacement of the oscillating body as a function of time.

Symbol:
    :code:`x(t)`
"""

time = Symbol("time", units.time, positive=True)
"""
Time.

Symbol:
    :code:`t`
"""

undamped_angular_frequency = Symbol("undamped_angular_frequency",
    angle_type / units.time,
    positive=True)
r"""
Undamped angular frequency of the oscillator.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

damping_ratio = Symbol("damping_ratio", dimensionless, positive=True)
r"""
Damping ratio, which critically determines the behavior of the system.

Symbol:
    :code:`z`

Latex:
    :math:`\zeta`
"""

definition = (Derivative(displacement(time), time, 2) +
    2 * damping_ratio * undamped_angular_frequency * Derivative(displacement(time), time) +
    undamped_angular_frequency**2 * displacement(time))
r"""
:code:`Derivative(x(t), (t, 2)) + 2 * z * w * Derivative(x(t), t) + w^2 * x(t) = 0`

Latex:
    .. math::
        \frac{d^2 x}{dt^2} + 2 \zeta \omega \frac{d x}{d t} + \omega^2 x(t) = 0
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
