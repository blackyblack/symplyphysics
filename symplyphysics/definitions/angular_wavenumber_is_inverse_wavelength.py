"""
Angular wavenumber is inverse wavelength
========================================

*Angular wavenumber* is the spatial frequency of a wave, measured in radians per unit distance.
"""

from sympy import Eq, pi
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    angle_type,
)

angular_wavenumber = Symbol("angular_wave_number", angle_type / units.length)
"""
Angular wavenumber of the wave.

Symbol:
    :code:`k`
"""

wavelength = Symbol("wavelength", units.length)
r"""
Wavelength of the wave.

Symbol:
    :code:`lambda`

Latex:
    :math:`\lambda`
"""

definition = Eq(angular_wavenumber, 2 * pi / wavelength)
r"""
:code:`k = 2 * pi / lambda`

Latex:
    .. math::
        k = \frac{2 \pi}{\lambda}
"""


@validate_input(wavelength_=wavelength)
@validate_output(angular_wavenumber)
def calculate_wavenumber(wavelength_: Quantity) -> Quantity:
    result = definition.rhs.subs(wavelength, wavelength_)
    return Quantity(result)
