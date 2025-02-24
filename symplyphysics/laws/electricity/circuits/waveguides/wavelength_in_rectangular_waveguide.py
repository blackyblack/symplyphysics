"""
Wavelength in rectangular waveguide
===================================

Wavelength in a rectangular waveguide is a function of operating wavelength (or
wavelength in free space) and the critical (or lower cutoff) wavelength, and is
always longer than that in free space.

**Links:**

#. `Microwaves101, section "Guide Wavelength" <https://www.microwaves101.com/encyclopedias/waveguide-mathematics>`__.
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

waveguide_wavelength = clone_as_symbol(symbols.wavelength, display_symbol="lambda_g", display_latex="\\lambda_\\text{g}")
"""
Guide :symbols:`wavelength` is defined as the distance between the two equal phase planes along the
waveguide.
"""

vacuum_wavelength = symbols.wavelength
"""
:symbols:`wavelength` of the signal wave in vacuum.
"""

critical_wavelength = clone_as_symbol(symbols.wavelength, display_symbol="lambda_c", display_latex="\\lambda_\\text{c}")
"""
Critical :symbols:`wavelength` of the waveguide. See :ref:`Critical wavelength of waveguide <critical_wavelength_waveguide_def>`.
"""

law = Eq(waveguide_wavelength,
    vacuum_wavelength / sqrt(1 - (vacuum_wavelength / critical_wavelength)**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(signal_wavelength_=vacuum_wavelength, critical_wavelength_=critical_wavelength)
@validate_output(waveguide_wavelength)
def calculate_waveguide_wavelength(signal_wavelength_: Quantity,
    critical_wavelength_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, waveguide_wavelength, dict=True)[0][waveguide_wavelength]
    result_expr = result_velocity_expr.subs({
        vacuum_wavelength: signal_wavelength_,
        critical_wavelength: critical_wavelength_
    })
    return Quantity(result_expr)
