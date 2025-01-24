"""
Wavelength in rectangular waveguide
===================================

A rectangular waveguide is a rectangular metal waveguide capable of supporting waves
propagating along it. There is a critical wavelength. Signals with a wavelength greater
than the critical one are attenuated and do not propagate in the waveguide.

"""

from sympy import Eq, solve, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

waveguide_wavelength = clone_as_symbol(symbols.wavelength, subscript="\\text{w}")
"""
:symbols:`wavelength` in the waveguide.
"""

signal_wavelength = symbols.wavelength
"""
:symbols:`wavelength` of the signal.
"""

critical_wavelength = clone_as_symbol(symbols.wavelength, subscript="\\text{c}")
"""
Critical :symbols:`wavelength` of the waveguide.
"""

law = Eq(waveguide_wavelength,
    signal_wavelength / sqrt(1 - (signal_wavelength / critical_wavelength)**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(signal_wavelength_=signal_wavelength, critical_wavelength_=critical_wavelength)
@validate_output(waveguide_wavelength)
def calculate_waveguide_wavelength(signal_wavelength_: Quantity,
    critical_wavelength_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, waveguide_wavelength, dict=True)[0][waveguide_wavelength]
    result_expr = result_velocity_expr.subs({
        signal_wavelength: signal_wavelength_,
        critical_wavelength: critical_wavelength_
    })
    return Quantity(result_expr)
