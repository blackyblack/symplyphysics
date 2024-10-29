"""
Peak wavelength via temperature
===============================

The wavelength of a blackbody's radiation peak is inversely proportional to its
temperature. This law is known as the *Wien's displacement law*.
"""

from sympy import Eq, solve
from symplyphysics import (
    symbols,
    quantities,
    Quantity,
    validate_input,
    validate_output,
    clone_as_symbol,
)

peak_wavelength = clone_as_symbol(
    symbols.wavelength,
    display_symbol="lambda_peak",
    display_latex="\\lambda_\\text{peak}",
)
"""
:symbols:`wavelength` at which the peak of the spectral radiance of a blackbody occurs.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the blackbody.
"""

law = Eq(peak_wavelength, quantities.wien_displacement_constant / temperature)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(object_temperature_=temperature)
@validate_output(peak_wavelength)
def calculate_intensive_wavelength(object_temperature_: Quantity) -> Quantity:
    result_wavelength_expr = solve(law, peak_wavelength, dict=True)[0][peak_wavelength]
    result_expr = result_wavelength_expr.subs({temperature: object_temperature_})
    return Quantity(result_expr)
