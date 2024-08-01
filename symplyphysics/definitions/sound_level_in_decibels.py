"""
Sound level in decibels
=======================

The sound level of a sound wave is a physical quantity that is used to describe the wave's
intensity.
"""

from sympy import Eq, log
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

sound_level = Symbol("sound_level", dimensionless)
r"""
Sound level of the sound wave.

Symbol:
    :code:`beta`

Latex:
    :math:`\beta`
"""

intensity = Symbol("intensity", units.power / units.area)
"""
Intensity of the sound wave.

Symbol:
    :code:`I`
"""

reference_sound_level = Quantity(10)
r"""
The sound level when the wave's intensity equals the reference intensity.

Symbol:
    :code:`beta0`

Latex:
    :math:`\beta_0`
"""

reference_intensity = Quantity(1e-12 * units.watt / units.meter**2)
"""
The intensity of a sound wave, relative to which the sound level is measured.

Symbol:
    :code:`I0`

Latex:
    :math:`I_0`
"""

definition = Eq(sound_level, reference_sound_level * log(intensity / reference_intensity, 10))
r"""
:code:`beta = beta0 * log10(I / I0)`

Latex:
    .. math::
        \beta = \beta_0 \log_{10} \left( \frac{I}{I_0} \right)
"""


@validate_input(intensity_=intensity)
@validate_output(sound_level)
def calculate_sound_level(intensity_: Quantity) -> float:
    result = definition.rhs.subs(intensity, intensity_)
    return Quantity(result).scale_factor
