"""
Angular wavenumber is inverse wavelength
========================================

*Angular wavenumber* is the spatial frequency of a wave, measured in radians per unit distance.
"""

from sympy import Eq, pi
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    angle_type,
)

angular_wavenumber = SymbolNew("k", angle_type / units.length)
"""
Angular wavenumber of the wave.
"""

wavelength = SymbolNew("lambda", units.length, display_latex="\\lambda")
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
