"""
Sound level in decibels
=======================

The sound level of a sound wave is a physical quantity that is used to describe the wave's
intensity.
"""

from sympy import Eq, log
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

sound_level = symbols.sound_intensity_level
"""
:symbols:`sound_intensity_level` of the sound wave.
"""

intensity = symbols.intensity
"""
:symbols:`intensity` of the sound wave.
"""

reference_sound_level = Quantity(10, display_symbol="L_I0", display_latex="L_{I0}")
"""
The sound level when the wave's intensity equals the reference intensity.
"""

reference_intensity = Quantity(1e-12 * units.watt / units.meter**2, display_symbol="I_0")
"""
The intensity of a sound wave relative to which the sound level is measured.
"""

definition = Eq(sound_level, reference_sound_level * log(intensity / reference_intensity, 10))
r"""
.. only:: comment

    Manual formulas are required since sympy does not handle logarithms with arbitrary bases well.

:code:`beta = beta0 * log_10(I / I0)`

Latex:
    .. math::
        \beta = \beta_0 \log_{10} \left( \frac{I}{I_0} \right)
"""


@validate_input(intensity_=intensity)
@validate_output(sound_level)
def calculate_sound_level(intensity_: Quantity) -> float:
    result = definition.rhs.subs(intensity, intensity_)
    return Quantity(result).scale_factor
