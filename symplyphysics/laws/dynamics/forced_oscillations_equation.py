"""
Forced oscillations equation
============================

*Forced, or driven, oscillations* are a type of oscillations in the presence of an external driving
force acting on the oscillating system. In the case of an oscillating external force, two angular
frequencies are associated with such a system: (1) the natural angular frequency of the system,
which is the angular frequency the system would oscillate with if no external force were present,
and (2) the angular frequency of the external force driving the oscillations. Such systems can
undergo resonance if the angular frequency of the driving force is close to the natural angular
frequency of the oscillator.
"""

from sympy import Derivative, Eq, cos, dsolve
from symplyphysics import (
    symbols,
    units,
    angle_type,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
)
from symplyphysics.core.symbols.quantities import scale_factor

displacement = Function("displacement", units.length)
"""
The displacement of the oscillating body from rest value.

Symbol:
    :code:`q(t)`
"""

mass = symbols.mass
"""
The :symbols:`mass` of the oscillating body.
"""

natural_angular_frequency = Symbol("natural_angular_frequency", angle_type / units.time)
r"""
The natural angular frequency of the oscillator.

Symbol:
    :code:`w0`

Latex:
    :math:`\omega_0`
"""

driving_force_amplitude = symbols.force
"""
The amplitude of the driving :symbols:`force`.
"""

driving_angular_frequency = Symbol("driving_angular_frequency", angle_type / units.time)
r"""
The angular frequency of the driving force.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

driving_phase_lag = Symbol("driving_phase_lag", angle_type)
r"""
The phase lag of the driving force.

Symbol:
    :code:`phi`

Latex:
    :math:`\varphi`
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

law = Eq(
    Derivative(displacement(time), time, 2) + natural_angular_frequency**2 * displacement(time),
    (driving_force_amplitude / mass) * cos(driving_angular_frequency * time + driving_phase_lag),
)
r"""
Derivative(q(t), (t, 2)) + w0^2 * q(t) = (F / m) * cos(w * t + phi)

Latex:
    .. math::
        \frac{d^2}{d t^2} q(t) + \omega_0^2 q(t) = \frac{F}{m} \cos{\left( \omega t + \varphi \right)}
"""


#pylint: disable=too-many-arguments
@validate_input(
    initial_position_=displacement,
    initial_velocity_=units.velocity,
    oscillator_mass_=mass,
    natural_angular_frequency_=natural_angular_frequency,
    driving_force_amplitude_=driving_force_amplitude,
    driving_angular_frequency_=driving_angular_frequency,
    driving_phase_lag_=driving_phase_lag,
    time_=time,
)
@validate_output(displacement)
def calculate_displacement(
    initial_position_: Quantity,
    initial_velocity_: Quantity,
    oscillator_mass_: Quantity,
    natural_angular_frequency_: Quantity,
    driving_force_amplitude_: Quantity,
    driving_angular_frequency_: Quantity,
    driving_phase_lag_: Quantity | float,
    time_: Quantity,
) -> Quantity:
    initial_conditions = {
        displacement(0): initial_position_,
        displacement(time).diff(time).subs(time, 0): initial_velocity_,
    }
    eqn = law.subs({
        natural_angular_frequency: natural_angular_frequency_,
        driving_angular_frequency: driving_angular_frequency_,
        mass: oscillator_mass_,
        driving_force_amplitude: driving_force_amplitude_,
        driving_phase_lag: scale_factor(driving_phase_lag_),
    })
    dsolved = dsolve(
        eqn,
        displacement(time),
        ics=initial_conditions,
    ).rhs
    result = dsolved.subs(time, time_)
    return Quantity(result)
