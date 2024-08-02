r"""
Wave equation in one dimension
==============================

The *wave equation* is a second-order linear partial differential equation used to
describe the propagation of waves, including standing wave fields such as mechanical
or electromagnetic waves.

**Notes:**

#. This equation is called one-dimensional because the displacement function depends
   only on one spatial dimension.
"""

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
    Quantity,
    validate_input,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.quantity_decorator import validate_output_same

displacement = symbols("displacement", cls=SymFunction, real=True)
"""
Factor representing a displacement from rest position, which could be
pressure, position, electric field, etc., as a function of position
and time.

Symbol:
    :code:`u(x, t)`
"""

position = Symbol("position", units.length, real=True)
"""
Position, or spatial variable.

Symbol:
    :code:`x`
"""

time = Symbol("time", units.time, positive=True)
"""
Time.

Symbol:
    :code:`t`
"""

phase_speed = Symbol("phase_speed", units.velocity, real=True)
"""
:doc:`Phase speed <laws.waves.phase_velocity_from_angular_frequency_and_wavenumber>` of the wave.

Symbol:
    :code:`v`
"""

definition = Eq(
    Derivative(displacement(position, time), position, 2),
    Derivative(displacement(position, time), time, 2) / phase_speed**2,
)
r"""
:code:`Derivative(u(x, t), (x, 2)) = (1/v^2) * Derivative(u(x, t), (t, 2))`

Latex:
    .. math::
        \frac{\partial^2 u}{\partial x^2} = \frac{1}{v^2} \frac{\partial^2 u}{\partial t^2}
"""

# A subset of solutions has the form `u = u_m * cos(x + v*t + phi)`
## u_m - amplitude of displacement
## phi - phase lag

amplitude = symbols("amplitude", positive=True)
phase_lag = symbols("phase_lag", real=True)
_length_unit = Symbol("length_unit", units.length)

# `phase_speed` can be negative or positive, depending on the direction of wave propagation
# - negative values denote propagation in the positive direction of x-axis
# - positive values denote propagation in the negative direction of x-axis
solution = amplitude * cos((position + phase_speed * time) / _length_unit + phase_lag)

_lhs = definition.lhs.subs(displacement(position, time), solution).doit()
_rhs = definition.rhs.subs(displacement(position, time), solution).doit()
assert expr_equals(_lhs, _rhs)

# Check that the equation holds for waves traveling in the opposite direction
_solution = solution.subs(phase_speed, -1 * phase_speed)
_lhs = definition.lhs.subs(displacement(position, time), _solution).doit()
_rhs = definition.rhs.subs(displacement(position, time), _solution).doit()
assert expr_equals(_lhs, _rhs)


@validate_input(
    phase_speed_=phase_speed,
    position_=position,
    time_=time,
)
@validate_output_same("amplitude_")
def calculate_displacement(
    amplitude_: Quantity,
    phase_speed_: Quantity,
    phase_lag_: float,
    position_: Quantity,
    time_: Quantity,
) -> Quantity:
    # This function only works for a subset of possible solutions of the wave equation,
    # namely those proportional to `cos(x + v*t + phi)` for some phase lag `phi`.

    result = solution.subs({
        amplitude: amplitude_,
        phase_speed: phase_speed_,
        phase_lag: phase_lag_,
        position: position_,
        time: time_,
        _length_unit: Quantity(1.0 * units.meter),
    })
    return Quantity(result)
