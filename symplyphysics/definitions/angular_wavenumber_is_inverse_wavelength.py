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

angular_wavenumber = Symbol("angular_wave_number", angle_type / units.length, display_symbol="k")
"""
Angular wavenumber of the wave.

Symbol:
    :code:`k`
"""

wavelength = Symbol("wavelength", units.length, display_symbol="lambda", display_latex="\\lambda")
"""
Wavelength of the wave.
"""

definition = Eq(angular_wavenumber, 2 * pi / wavelength)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(wavelength_=wavelength)
@validate_output(angular_wavenumber)
def calculate_wavenumber(wavelength_: Quantity) -> Quantity:
    result = definition.rhs.subs(wavelength, wavelength_)
    return Quantity(result)
