r"""
Displacement in forced non-resonant oscillations
================================================

*Forced, or driven, oscillations* are a type of oscillations in the presence of an external driving
force acting on the oscillating system. In the case of an oscillating external force, two angular
frequencies are associated with such a system: (1) the natural angular frequency of the system,
which is the angular frequency the system would oscillate with if no external force were present,
and (2) the angular frequency of the external force driving the oscillations.

**Conditions:**

#. Angular frequency of the external force is strictly not equal to the natural angular frequency of the oscillator.
#. No damping is present in the system.

**Notes:**

#. The external driving force has the form of :math:`f(t) = f_m \cos{\left( \omega t + \varphi \right)}`.
#. The complete expression of the displacement function can be found as the sum of the solution of 
   simple harmonic motion equation and the particular solution presented here.

**Links:**

#. `Physics LibreTexts, derivable from (15.7.2) and (15.7.3) <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/Book%3A_University_Physics_I_-_Mechanics_Sound_Oscillations_and_Waves_(OpenStax)/15%3A_Oscillations/15.07%3A_Forced_Oscillations>`__.
"""

from sympy import Eq, cos, dsolve
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
    clone_as_symbol,
    clone_as_function,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.symbols.quantities import scale_factor
from symplyphysics.laws.dynamics import forced_oscillations_equation as forced_eqn

time = symbols.time
"""
:symbols:`time`.
"""

displacement = clone_as_function(symbols.position, [time], display_symbol="q", display_latex="q")
"""
The particular solution of the forced oscillations equation that accounts for the
oscillator's response to the driving force. See :symbols:`position`.
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
"""
The natural :symbols:`angular_frequency` of the oscillator.
"""

driving_force_amplitude = symbols.force
"""
The amplitude of the external driving :symbols:`force`.
"""

driving_angular_frequency = symbols.angular_frequency
"""
The :symbols:`angular_frequency` of the external driving force.
"""

driving_phase_lag = symbols.phase_shift
"""
The :symbols:`phase_shift` of the oscillations of the external force.
"""

law = Eq(displacement(time),
    (driving_force_amplitude / mass) * cos(driving_angular_frequency * time + driving_phase_lag) /
    (natural_angular_frequency**2 - driving_angular_frequency**2))
"""
:laws:symbol::

:laws:latex::
"""

# Derive law from forced oscillations equation

_eqn = forced_eqn.law.subs({
    forced_eqn.mass: mass,
    forced_eqn.natural_angular_frequency: natural_angular_frequency,
    forced_eqn.driving_force_amplitude: driving_force_amplitude,
    forced_eqn.driving_angular_frequency: driving_angular_frequency,
    forced_eqn.driving_phase_lag: driving_phase_lag,
    forced_eqn.time: time,
})

_dsolved = dsolve(_eqn, forced_eqn.displacement(time)).rhs

# We're interested in the particular solution of the differential equation, therefore
# throwing out the part that describes the oscillator's motion in the absence of an
# external driving force.
_dsolved = _dsolved.subs({"C1": 0, "C2": 0})

assert expr_equals(_dsolved, law.rhs)


@validate_input(
    oscillator_mass_=mass,
    natural_angular_frequency_=natural_angular_frequency,
    driving_force_amplitude_=driving_force_amplitude,
    driving_angular_frequency_=driving_angular_frequency,
    driving_phase_lag_=driving_phase_lag,
    time_=time,
)
@validate_output(displacement)
def calculate_driven_displacement(
    oscillator_mass_: Quantity,
    natural_angular_frequency_: Quantity,
    driving_force_amplitude_: Quantity,
    driving_angular_frequency_: Quantity,
    driving_phase_lag_: Quantity | float,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        mass: oscillator_mass_,
        natural_angular_frequency: natural_angular_frequency_,
        driving_force_amplitude: driving_force_amplitude_,
        driving_angular_frequency: driving_angular_frequency_,
        driving_phase_lag: scale_factor(driving_phase_lag_),
        time: time_,
    })
    return Quantity(result)
