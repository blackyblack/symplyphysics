"""
Angular wavenumber is inverse wavelength
========================================

*Angular wavenumber* is the spatial frequency of a wave, measured in radians per unit distance.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Wavenumber#Definition>`__.
"""

from sympy import Eq, pi
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

angular_wavenumber = symbols.angular_wavenumber
"""
:symbols:`angular_wavenumber` of the wave.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength` of the wave.
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
