"""
Characteristic resitance of rectangular waveguide for transverse magnetic waves
===============================================================================

The resistance of a rectangular waveguide to transverse magnetic waves can be calculated
from the resistance of the medium within the waveguide, the wavelength of the signal
and the critical wavelength of the waveguide.

**Conditions:**

#. Waves propagating in the waveguide must be transverse magnetic waves.

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

wave_impedance = symbols.wave_impedance
"""
:symbols:`wave_impedance` in the waveguide.
"""

medium_impedance = clone_as_symbol(symbols.wave_impedance, subscript="0")
"""
:symbols:`wave_impedance` of the medium filling the waveguide.
"""

vacuum_wavelength = symbols.wavelength
"""
:symbols:`wavelength` of the signal in vacuum.
"""

critical_wavelength = clone_as_symbol(symbols.wavelength, display_symbol="lambda_c", display_latex="\\lambda_\\text{c}")
"""
Critical :symbols:`wavelength`. See :ref:`Critical wavelength of waveguide`.
"""

law = Eq(wave_impedance, medium_impedance * sqrt(1 - (vacuum_wavelength / critical_wavelength)**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(resistance_of_medium_=wave_impedance,
    wavelength_=vacuum_wavelength,
    critical_wavelength_=critical_wavelength)
@validate_output(wave_impedance)
def calculate_resistance(resistance_of_medium_: Quantity, wavelength_: Quantity,
    critical_wavelength_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, wave_impedance, dict=True)[0][wave_impedance]
    result_expr = result_velocity_expr.subs({
        medium_impedance: resistance_of_medium_,
        vacuum_wavelength: wavelength_,
        critical_wavelength: critical_wavelength_
    })
    return Quantity(result_expr)
