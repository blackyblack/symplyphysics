r"""
Displacement in resonant oscillations
=====================================

When an oscillating external force is driving the oscillations of an oscillator,
amplitude of oscillations is greatest when the angular frequency of the driving
is equal to the natural frequency of the oscillator. This condition is called
resonance.

**Conditions:**

#. Angular frequency of the driving force is equal to the natural frequency of the oscillator.
#. No damping is present in the system.

**Notes:**

#. The expression of the driving force has the form :math:`F \cos{\left( \omega t + \varphi \right)}`
   where :math:`\omega` is the angular frequency of its oscillations.
"""

from sympy import Eq, dsolve, sin
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.symbols.quantities import scale_factor
from symplyphysics.laws.dynamics import forced_oscillations_equation as forced_eqn

time = symbols.time
"""
:symbols:`time`.
"""

resonant_displacement = clone_as_function(symbols.position, [time])
"""
The displacement of resonant oscillations as a function of :attr:`~time`. See :symbols:`position`.
"""

mass = symbols.mass
"""
The :symbols:`mass` of the oscillator.
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
The amplitude of the driving :symbols:`force`.
"""

driving_phase_lag = symbols.phase_shift
"""
The :symbols:`phase_shift` of the oscillations of the driving force.
"""

law = Eq(resonant_displacement(time), (driving_force_amplitude / mass) * time *
    sin(natural_angular_frequency * time + driving_phase_lag) / (2 * natural_angular_frequency))
"""
:laws:symbol::

:laws:latex::
"""

# Derive law from driven oscillations equation

_eqn = forced_eqn.law.subs({
    forced_eqn.mass: mass,
    forced_eqn.natural_angular_frequency: natural_angular_frequency,
    forced_eqn.driving_force_amplitude: driving_force_amplitude,
    forced_eqn.driving_angular_frequency: natural_angular_frequency,
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
    driving_phase_lag_=driving_phase_lag,
    time_=time,
)
@validate_output(resonant_displacement)
def calculate_resonant_displacement(
    oscillator_mass_: Quantity,
    natural_angular_frequency_: Quantity,
    driving_force_amplitude_: Quantity,
    driving_phase_lag_: Quantity | float,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        mass: oscillator_mass_,
        natural_angular_frequency: natural_angular_frequency_,
        driving_force_amplitude: driving_force_amplitude_,
        driving_phase_lag: scale_factor(driving_phase_lag_),
        time: time_,
    })
    return Quantity(result)
