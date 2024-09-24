"""
Forced oscillations equation
============================

*Forced, or driven, oscillations* are a type of oscillations in the presence of an external driving
force acting on the oscillating system. In the case of an oscillating external force, two angular
frequencies are associated with such a system: 

#. the *natural angular frequency* of the system, which is the angular frequency the system would 
   oscillate with if no external force were present,

#. the angular frequency of the external force driving the oscillations.

Such systems can undergo resonance if the angular frequency of the driving force is close to the
natural angular frequency of the oscillator.
"""

from sympy import Derivative, Eq, cos, dsolve
from symplyphysics import (
    symbols,
    units,
    Quantity,
    validate_input,
    validate_output,
    clone_as_function,
    clone_as_symbol,
)
from symplyphysics.core.symbols.quantities import scale_factor

displacement = clone_as_function(symbols.position, display_symbol="x(t)")
"""
The displacement of the oscillating body from rest value. See :symbols:`position`.
"""

mass = symbols.mass
"""
The :symbols:`mass` of the oscillating body.
"""

natural_angular_frequency = clone_as_symbol(
    symbols.angular_frequency,
    display_symbol="w_0",
    display_latex="\\omega_0",
)
r"""
The natural :symbols:`angular_frequency` of the oscillator.
"""

driving_force_amplitude = symbols.force
"""
The amplitude of the driving :symbols:`force`.
"""

driving_angular_frequency = symbols.angular_frequency
r"""
The :symbols:`angular_frequency` of the driving force.
"""

driving_phase_lag = symbols.phase_shift
r"""
The :symbols:`phase_shift` of the driving force.
"""

time = symbols.time
"""
:symbols:`time`.
"""

law = Eq(
    Derivative(displacement(time), time, 2) + natural_angular_frequency**2 * displacement(time),
    (driving_force_amplitude / mass) * cos(driving_angular_frequency * time + driving_phase_lag),
)
"""
:laws:symbol::

:laws:latex::
"""


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
    # pylint: disable=too-many-arguments, too-many-positional-arguments
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
