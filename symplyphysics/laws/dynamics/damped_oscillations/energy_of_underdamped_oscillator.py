"""
Energy of underdamped oscillator
================================

In the presence of a damping force the oscillating system is no longer closed and its
energy dissipates to the environment. The total energy of the oscillator becomes
converted into thermal energy. For small values of the damping ratio, the equation given
in this law approximately describes the total mechanical energy of the underdamped
oscillator.

**Conditions:**

#. Damping ratio :math:`\\zeta` is small, i.e. :math:`\\zeta \\ll 1`.

**Links:**

#. Similar equation (15-44) on p. 431 of "Fundamentals of Physics" by David Halladay et al., 10th Ed.
"""

from sympy import Eq, exp
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    clone_as_symbol,
)

energy = symbols.energy
"""
:symbols:`energy` of the oscillator.
"""

mass = symbols.mass
"""
:symbols:`mass` of the oscillator.
"""

amplitude = Symbol("A", units.length)
"""
Amplitude, or maximum absolute displacement, of the oscillator.
"""

undamped_angular_frequency = clone_as_symbol(symbols.angular_frequency, subscript="0")
"""
:symbols:`angular_frequency` of the undamped oscillator.
"""

exponential_decay_constant = symbols.exponential_decay_constant
"""
:symbols:`exponential_decay_constant` of the oscilaltor.
"""

time = symbols.time
"""
:symbols:`time`.
"""

law = Eq(
    energy, mass * undamped_angular_frequency**2 * amplitude**2 *
    exp(-2 * exponential_decay_constant * time) / 2)
"""
:laws:symbol::

:laws:latex::
"""

# TODO Derive from [underdamped oscillations](../../kinematics/damped_oscillations/underdamping.py)


@validate_input(
    mass_=mass,
    maximum_amplitude_=amplitude,
    undamped_angular_frequency_=undamped_angular_frequency,
    exponential_decay_constant_=exponential_decay_constant,
    time_=time,
)
@validate_output(energy)
def calculate_oscillator_energy(
    mass_: Quantity,
    maximum_amplitude_: Quantity,
    undamped_angular_frequency_: Quantity,
    exponential_decay_constant_: Quantity,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        mass: mass_,
        amplitude: maximum_amplitude_,
        undamped_angular_frequency: undamped_angular_frequency_,
        exponential_decay_constant: exponential_decay_constant_,
        time: time_,
    })
    return Quantity(result)
