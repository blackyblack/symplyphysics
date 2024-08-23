r"""
Displacement in critical damping
================================

In the presence of a damping force in the oscillating system, the system's behaviour
depends on the value of the damping ratio. When it is equal to :math:`1`, the system returns
to equilibrium as quickly as possible without oscillating, but overshoot can occur if
initial velocity is nonzero. This behaviour is also called *critical damping*.

**Conditions:**

#. The system is critically damped, i.e. its damping ratio :math:`\zeta = 1`.
"""

from sympy import Eq, exp, dsolve, solve, Function as SymFunction, symbols
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    angle_type,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import damped_harmonic_oscillator_equation as damped_eqn

displacement = Symbol("displacement", units.length, real=True)
"""
Displacement from rest, usually a function of time.

Symbol:
    :code:`q`
"""

time = Symbol("time", units.time, nonnegative=True)
r"""
Time at which :math:`q` is measured.

Symbol:
    :code:`t`
"""

undamped_angular_frequency = Symbol("undamped_angular_frequency",
    angle_type / units.time,
    positive=True)
r"""
Angular frequency of the undamped oscillator.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

initial_position = Symbol("initial_position", units.length, real=True)
r"""
Initial position of the oscillator.

Symbol:
    :code:`x_0`

Latex:
    :math:`x_0`
"""

initial_speed = Symbol("initial_speed", units.velocity, real=True)
r"""
Initial speed of the oscillator.

Symbol:
    :code:`v_0`

Latex:
    :math:`v_0`
"""

law = Eq(
    displacement,
    exp(-1 * undamped_angular_frequency * time) * (initial_position +
    (initial_speed + initial_position * undamped_angular_frequency) * time),
)
r"""
:code:`q = exp(-1 * w * t) * (x_0 + (v_0 + x_0 * w) * t)`

Latex:
    .. math::
        q = \exp(- \omega t) (x_0 + (v_0 + x_0 \omega) t)
"""

# Derive from damped oscillator equation

_displacement = symbols("displacement", real=True, cls=SymFunction)

_eqn = damped_eqn.definition.subs(damped_eqn.time, time).subs({
    damped_eqn.displacement(time): _displacement(time),
    damped_eqn.undamped_angular_frequency: undamped_angular_frequency,
    damped_eqn.damping_ratio: 1,
})

_dsolved = dsolve(_eqn, _displacement(time)).rhs

_initial_position_eqn = Eq(initial_position, _dsolved.subs(time, 0))
_initial_velocity_eqn = Eq(initial_speed, _dsolved.diff(time).subs(time, 0))

_c12 = solve(
    [_initial_position_eqn, _initial_velocity_eqn],
    ("C1", "C2"),
    dict=True,
)[0]

_dsolved = _dsolved.subs(_c12)

assert expr_equals(law.rhs, _dsolved)


@validate_input(
    initial_position_=initial_position,
    initial_velocity_=initial_speed,
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
        initial_speed: initial_velocity_,
        undamped_angular_frequency: undamped_angular_frequency_,
        time: time_,
    })
    return Quantity(result)
