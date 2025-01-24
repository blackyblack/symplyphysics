"""
Characteristic resitance of rectangular waveguide for transverse electric waves
===============================================================================

A rectangular waveguide is a rectangular metal waveguide capable of supporting waves
propagating along it. There is a critical wavelength. Signals with a wavelength greater
than the critical one are attenuated and do not propagate in the waveguide.

**Conditions:**

#. Waves propagating in the waveguide must be transverse electric waves.

..
    TODO: find link
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

resistance = symbols.electrical_resistance
"""
:symbols:`electrical_resistance` of the waveguide.
"""

medium_resistance = clone_as_symbol(symbols.electrical_resistance, subscript="0")
"""
:symbols:`electrical_resistance` of the material filling the waveguide.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength`.
"""

critical_wavelength = clone_as_symbol(symbols.wavelength, subscript="\\text{c}")
"""
Critical :symbols:`wavelength`.
"""

law = Eq(resistance, medium_resistance / sqrt(1 - (wavelength / critical_wavelength)**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(resistance_of_medium_=medium_resistance,
    wavelength_=wavelength,
    critical_wavelength_=critical_wavelength)
@validate_output(resistance)
def calculate_resistance(resistance_of_medium_: Quantity, wavelength_: Quantity,
    critical_wavelength_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, resistance, dict=True)[0][resistance]
    result_expr = result_velocity_expr.subs({
        medium_resistance: resistance_of_medium_,
        wavelength: wavelength_,
        critical_wavelength: critical_wavelength_
    })
    return Quantity(result_expr)
