"""
Section length of multistage transformer
========================================

A multistage resistance transformer consists of several sections. The length of each
section depends on the wavelength at the upper operating frequency and the wavelength at
the lower operating frequency.

..
    TODO: find link
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

section_length = symbols.length
"""
:symbols:`length` of the section.
"""

upper_frequency_wavelength = clone_as_symbol(symbols.wavelength, subscript="1")
"""
:symbols:`wavelength` at the upper operating :symbols:`temporal_frequency`.
"""

lower_frequency_wavelength = clone_as_symbol(symbols.wavelength, subscript="2")
"""
:symbols:`wavelength` at the lower operating :symbols:`temporal_frequency`.
"""

law = Eq(
    section_length, upper_frequency_wavelength * lower_frequency_wavelength / (2 *
    (upper_frequency_wavelength + lower_frequency_wavelength)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(wavelength_for_upper_frequency_=upper_frequency_wavelength,
    lower_frequency_wavelength_=lower_frequency_wavelength)
@validate_output(section_length)
def calculate_length_of_section(wavelength_for_upper_frequency_: Quantity,
    lower_frequency_wavelength_: Quantity) -> Quantity:
    result_expr = solve(law, section_length, dict=True)[0][section_length]
    result_expr = result_expr.subs({
        upper_frequency_wavelength: wavelength_for_upper_frequency_,
        lower_frequency_wavelength: lower_frequency_wavelength_,
    })
    return Quantity(result_expr)
