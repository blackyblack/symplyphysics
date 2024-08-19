r"""
Resonant frequencies of stretched string with fixed ends
========================================================

For a string with fixed ends there is only a limited set of frequencies at which standing waves will
occur on it. Each possible frequency is a resonant frequency, and the corresponding wave pattern is
an oscillation mode. The oscillation mode corresponding to :math:`m = 1` is called the fundamental mode or
the first harmonic, the mode corresponding to :math:`m = 2` is the second harmonic, and so on.
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    dimensionless,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import (
    temporal_frequency_from_period as frequency_law,)
from symplyphysics.laws.waves import (
    wavelength_from_phase_speed_and_period as wavelength_law,
    wavelength_of_standing_wave_in_string_with_fixed_ends as resonance_law,
)

resonant_frequency = Symbol("resonant_frequency", units.frequency, positive=True)
r"""
Resonant frequency of the :math:`m`-th harmonic.

Symbol:
    :code:`f`
"""

harmonic_number = Symbol("harmonic_number", dimensionless, integer=True, positive=True)
"""
An integer called harmonic number.

Symbol:
    :code:`m`
"""

phase_velocity = Symbol("phase_velocity", units.velocity, positive=True)
"""
Phase velocity of the wave.

Symbol:
    :code:`v`
"""

string_length = Symbol("string_length", units.length, positive=True)
"""
Length of the string.

Symbol: 
    :code:`l`
"""

law = Eq(resonant_frequency, harmonic_number * phase_velocity / (2 * string_length))
r"""
:code:`f = m * v / (2 * l)`

Latex:
    .. math::
        f = \frac{m v}{2 l}
"""

# Derive from the condition for a standing wave along a stretched string with fixed ends

_wavelength = solve(resonance_law.law, resonance_law.wavelength)[0].subs({
    resonance_law.integer_factor: harmonic_number,
    resonance_law.string_length: string_length,
})

_period = solve(wavelength_law.law, wavelength_law.period)[0].subs({
    wavelength_law.wavelength: _wavelength,
    wavelength_law.phase_velocity: phase_velocity,
})

_frequency = solve(frequency_law.law,
    frequency_law.temporal_frequency)[0].subs(frequency_law.period, _period)

assert expr_equals(_frequency, law.rhs)


@validate_input(
    harmonic_number_=harmonic_number,
    phase_velocity_=phase_velocity,
    string_length_=string_length,
)
@validate_output(resonant_frequency)
def calculate_resonant_frequency(
    harmonic_number_: int,
    phase_velocity_: Quantity,
    string_length_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        harmonic_number: harmonic_number_,
        phase_velocity: phase_velocity_,
        string_length: string_length_,
    })
    return Quantity(result)
