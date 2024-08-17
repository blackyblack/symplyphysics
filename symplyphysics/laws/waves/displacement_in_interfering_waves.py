r"""
Displacement in interfering waves
=================================

If two waves are traveling in the same direction and have the same amplitude, period, and
wavelength (and hence the same frequency and wavenumber), but differ in the phase constant, the result
is a single wave with the same period and wavelength, but its amplitude depends on the phase
shift between the waves. If the shift is a multiple of :math:`2 \pi`, the waves are exactly
in phase and their interference is *fully constructive*. If it is :math:`\pi` plus a multiple of :math:`2 \pi`,
they are exactly out of phase and their interference is *fully destructive*.

**Notes:**

#. The form of the first wave is :math:`u_\text{max} \sin(k x - \omega t)` and of the second wave
   is :math:`u_\text{max} \sin(k x - \omega t + \varphi)`.
#. The travel of the waves in unaffected by their interference.

**Conditions:**

#. The waves are traveling in the same (or similar) directions.
#. They have the same amplitude, wavenumber and frequency.
"""

from sympy import Eq, sin, cos, symbols, Function as SymFunction
from symplyphysics import (
    units,
    angle_type,
    Symbol,
    Quantity,
    validate_input,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.symbols.quantities import scale_factor
from symplyphysics.core.quantity_decorator import validate_output_same

total_displacement = symbols("total_displacement")
"""
Displacement of the resulting wave.

Symbol:
    :code:`u`
"""

amplitude = symbols("amplitude")
r"""
Amplitude of the interfering waves.

Symbol:
    :code:`u_max`

Latex:
    :math:`u_\text{max}`
"""

phase_shift = Symbol("phase_shift", angle_type)
r"""
Phase shift between the interphering waves.

Symbol:
    :code:`phi`

Latex:
    :math:`\varphi`
"""

angular_wavenumber = Symbol("angular_wavenumber", angle_type / units.length)
r"""
Angular wavenumber of the interfering waves.

Symbol:
    :code:`k`
"""

angular_frequency = Symbol("angular_frequency", angle_type / units.time)
r"""
Angular frequency of the interfering waves.

Symbol:
    :code:`w`

Latex:
    :math:`\omega`
"""

position = Symbol("position", units.length)
"""
Position, or spatial coordinate.

Symbol:
    :code:`x`
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

law = Eq(total_displacement, (2 * amplitude) * cos(phase_shift / 2) *
    sin(angular_wavenumber * position - angular_frequency * time + phase_shift / 2))
r"""
:code:`u = 2 * u_max * cos(phi / 2) * sin(k * x - w * t + phi / 2)`

Latex:
    .. math::
        u = 2 u_\text{max} \cos \left( \frac{\varphi}{2} \right) \sin \left( k x - \omega t + \frac{\varphi}{2} \right)
"""

# Derive from the sum of two waves

# The form of the waves is taken from the notes above
_first_wave = amplitude * sin(angular_wavenumber * position - angular_frequency * time)
_second_wave = amplitude * sin(angular_wavenumber * position - angular_frequency * time +
    phase_shift)

_sum_of_waves = (_first_wave + _second_wave).simplify()

assert expr_equals(_sum_of_waves, law.rhs)


@validate_input(
    phase_shift_=phase_shift,
    angular_wavenumber_=angular_wavenumber,
    angular_frequency_=angular_frequency,
    position_=position,
    time_=time,
)
@validate_output_same("amplitude_")
def calculate_displacement(
    amplitude_: Quantity,
    phase_shift_: Quantity | float,
    angular_wavenumber_: Quantity,
    angular_frequency_: Quantity,
    position_: Quantity,
    time_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        amplitude: amplitude_,
        phase_shift: scale_factor(phase_shift_),
        angular_wavenumber: angular_wavenumber_,
        angular_frequency: angular_frequency_,
        position: position_,
        time: time_,
    })
    return Quantity(result)
