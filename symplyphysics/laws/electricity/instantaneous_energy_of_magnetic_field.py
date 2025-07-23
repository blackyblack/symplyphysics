"""
Instantaneous energy of magnetic field
======================================

There is an oscillatory circuit with alternating current. Then the energy of the magnetic field in
the inductor will depend on the inductance, the maximum value of the current, the angular frequency
of the current, the time and the initial phase.

**Conditions:**

#. The current depends on time sinusoidally:

   .. math::

       I(t) = I_\\text{max} \cos(\omega t + \\varphi)
"""

from sympy import Eq, solve, cos
from symplyphysics import Quantity, validate_input, validate_output, symbols, clone_as_symbol

energy = symbols.energy
"""
:symbols:`energy` stored in the coil.
"""

inductance = symbols.inductance
"""
:symbols:`inductance` of the coil.
"""

current_amplitude = clone_as_symbol(symbols.current,
    display_symbol="I_max",
    display_latex="I_\\text{max}")
"""
:symbols:`current` amplitude.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the current.
"""

time = symbols.time
"""
:symbols:`time`.
"""

initial_phase = symbols.phase_shift
"""
Initial :symbols:`phase_shift` of the oscillations.
"""

law = Eq(energy,
    (inductance * current_amplitude**2 / 2) * cos(angular_frequency * time + initial_phase)**2)
"""
:laws:symbol::

:laws:latex::
"""

# Derivable from the formula for energy stored in the inductor: W = L * I^2 / 2 (TODO: make this law)


@validate_input(inductance_=inductance,
    current_amplitude_=current_amplitude,
    frequency_=angular_frequency,
    time_=time,
    initial_phase_=initial_phase)
@validate_output(energy)
def calculate_energy(inductance_: Quantity, current_amplitude_: Quantity, frequency_: Quantity,
    time_: Quantity, initial_phase_: float | Quantity) -> Quantity:
    result_expr = solve(law, energy, dict=True)[0][energy]
    result_expr = result_expr.subs({
        inductance: inductance_,
        current_amplitude: current_amplitude_,
        angular_frequency: frequency_,
        time: time_,
        initial_phase: initial_phase_
    })
    return Quantity(result_expr)
