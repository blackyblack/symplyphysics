"""
Displacement in simple harmonic motion
======================================

Any motion that repeats at regular intervals is called periodic, or harmonic, motion.
Simple harmonic motion is a particular type of repeated motion that is a sinusoidal
function of time.
"""

from sympy import Eq, cos, symbols, dsolve, Function as SymFunction
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    angle_type,
    validate_input,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.quantity_decorator import validate_output_same
from symplyphysics.core.symbols.quantities import scale_factor
from symplyphysics.definitions import harmonic_oscillator_is_second_derivative_equation as harmonic_eqn

displacement = symbols("displacement", real=True)
"""
Displacement from rest, usually a function of time.

Symbol:
    :code:`q`
"""

amplitude = symbols("amplitude", positive=True)
r"""
The maximum absolute value of the displacement.

Symbol:
    :code:`q_max`

Latex:
    :math:`q_\text{max}`
"""

angular_frequency = Symbol("angular_frequency", angle_type / units.time, positive=True)
r"""
Angular frequency of oscillations.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

phase_shift = Symbol("phase_shift", angle_type, real=True)
r"""
Phase shift of oscillations, which is the phase at :math:`t = 0`.

Symbol:
    :code:`phi`

Latex:
    :math:`\varphi`
"""

time = Symbol("time", units.time, real=True)
"""
Time at which :math:`q` is measured.

Symbol:
    :code:`t`
"""

law = Eq(displacement, amplitude * cos(angular_frequency * time + phase_shift))
r"""
:code:`q = q_max * cos(w * t + phi)`

Latex:
    .. math::
        q = q_\text{max} \cos(\omega t + \varphi)
"""

# Derive law from [oscillator equation](../../definitions/harmonic_oscillator_is_second_derivative_equation.py)

_displacement = symbols("displacement", cls=SymFunction, real=True)

_eqn = harmonic_eqn.definition.replace(harmonic_eqn.displacement, _displacement).subs({
    harmonic_eqn.time: time,
    harmonic_eqn.angular_frequency: angular_frequency,
})

_initial_position = law.rhs.subs(time, 0)
_initial_velocity = law.rhs.diff(time).subs(time, 0)

_dsolved = dsolve(_eqn,
    _displacement(time),
    ics={
    _displacement(0): _initial_position,
    _displacement(time).diff(time).subs(time, 0): _initial_velocity,
    }).rhs

assert expr_equals(law.rhs, _dsolved)


@validate_input(
    angular_frequency_=angular_frequency,
    phase_lag_=phase_shift,
    time_=time,
)
@validate_output_same("amplitude_")
def calculate_displacement(
    amplitude_: Quantity,
    angular_frequency_: Quantity,
    phase_lag_: Quantity | float,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        amplitude: amplitude_,
        angular_frequency: angular_frequency_,
        phase_shift: scale_factor(phase_lag_),
        time: time_,
    })
    return Quantity(result)
