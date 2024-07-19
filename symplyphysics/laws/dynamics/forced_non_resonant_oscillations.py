r"""
Forced non-resonant oscillations
================================

*Forced, or driven, *oscillations* are a type of oscillations in the precence of an external driving
force acting on the oscillating system. In the case of an oscillating external force, two angular
frequencies are associated with such a system: (1) the natural angular frequency of the system,
which is the angular frequency the system would oscillate with if no external force were present,
and (2) the angular frequency of the external force driving the oscillations.

**Conditions:**

#. Angular frequency of the external force is strictly not equal to the natural angular frequency of the oscillator.
#. No damping is present in the system.

**Notes:**

#. The external driving force has the form of :math:`f(t) = f_m \cos{\omega t + \varphi}`.
#. The complete expression of the displacement function can be found as the sum of the solution of 
   simple harmonic motion equation and the particular solution presented here.
"""

from sympy import Eq, cos, dsolve
from symplyphysics import (
    clone_symbol,
    symbols,
    units,
    Quantity,
    Symbol,
    Function,
    angle_type,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.symbols.quantities import scale_factor
from symplyphysics.laws.dynamics import forced_oscillations_equation as forced_eqn

driven_displacement = Function("driven_displacement", units.length)
"""
The particular solution of the forced oscillations equation that accounts for the
oscillator's response to the driving force.

Symbol:
    q(t)
"""

oscillator_mass = clone_symbol(symbols.basic.mass, "oscillator_mass")
"""
The :attr:`~symplyphysics.symbols.basic.mass` of the oscillating body.

Symbol:
    m
"""

natural_angular_frequency = Symbol("natural_angular_frequency", angle_type / units.time)
r"""
The natural angular frequency of the oscillator.

Symbol:
    omega0

Latex:
    :math:`\omega_0`
"""

driving_force_amplitude = clone_symbol(symbols.dynamics.force, "driving_force_amplitude")
"""
The amplitude of the external driving :attr:`~symplyphysics.symbols.dynamics.force`.

Symbol:
    f
"""

driving_angular_frequency = Symbol("driving_angular_frequency", angle_type / units.time)
r"""
The angular frequency of the external driving force.

Symbol:
    omega

Latex:
    :math:`\omega`
"""

driving_phase_lag = Symbol("driving_phase_lag", angle_type)
r"""
The phase lag of the oscillations of the external force.

Symbol:
    phi

Latex:
    :math:`\varphi`
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    t
"""

law = Eq(driven_displacement(time), (driving_force_amplitude / oscillator_mass) *
    cos(driving_angular_frequency * time + driving_phase_lag) /
    (natural_angular_frequency**2 - driving_angular_frequency**2))
r"""
q(t) = (f / m) * cos(omega * t + phi) / (omega0**2 - omega**2)

Latex:
    :math:`q(t) = \frac{f}{m \left( \omega_0^2 - \omega^2 \right)} \cos{\left( \omega t + \varphi \right)}`
"""

# Derive law from forced oscillations equation

_eqn = forced_eqn.law.subs({
    forced_eqn.oscillator_mass: oscillator_mass,
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


#pylint: disable=too-many-arguments
@validate_input(
    oscillator_mass_=oscillator_mass,
    natural_angular_frequency_=natural_angular_frequency,
    driving_force_amplitude_=driving_force_amplitude,
    driving_angular_frequency_=driving_angular_frequency,
    driving_phase_lag_=driving_phase_lag,
    time_=time,
)
@validate_output(driven_displacement)
def calculate_driven_displacement(
    oscillator_mass_: Quantity,
    natural_angular_frequency_: Quantity,
    driving_force_amplitude_: Quantity,
    driving_angular_frequency_: Quantity,
    driving_phase_lag_: Quantity | float,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        oscillator_mass: oscillator_mass_,
        natural_angular_frequency: natural_angular_frequency_,
        driving_force_amplitude: driving_force_amplitude_,
        driving_angular_frequency: driving_angular_frequency_,
        driving_phase_lag: scale_factor(driving_phase_lag_),
        time: time_,
    })
    return Quantity(result)
