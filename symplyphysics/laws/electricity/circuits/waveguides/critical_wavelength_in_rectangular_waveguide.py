"""
Critical wavelength in rectangular waveguide
============================================

A rectangular waveguide is a rectangular metal waveguide capable of supporting waves
propagating along it. The critical frequency depends on the indices of the specific
propagation mode and the dimensions of the waveguide cross-section.

..
    TODO: find link
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    dimensionless,
    symbols,
    clone_as_symbol,
    Symbol,
)

critical_wavelength = clone_as_symbol(symbols.wavelength, display_symbol="lambda_c", display_latex="\\lambda_\\text{c}")
"""
Critical :symbols:`wavelength` in the waveguide. See :ref:`Critical wavelength of waveguide <critical_wavelength_waveguide_def>`.
"""

first_index = Symbol("m", dimensionless)
"""
The first index shows how many half-wavelengths fit across the width of the cross
section.
"""

second_index = Symbol("n", dimensionless)
"""
The second index shows how many half-wavelengths fit across the height of the cross
section.
"""

width = clone_as_symbol(symbols.length, display_symbol="a", display_latex="a")
"""
Width, or first dimension of the cross section. See :symbols:`length`.
"""

height = clone_as_symbol(symbols.length, display_symbol="b", display_latex="b")
"""
Height, or second dimension of the cross section. See :symbols:`length`.
"""

law = Eq(critical_wavelength, 2 / sqrt((first_index / width)**2 + (second_index / height)**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(first_index_=first_index, second_index_=second_index, width_=width, height_=height)
@validate_output(critical_wavelength)
def calculate_critical_wavelength(first_index_: float, second_index_: float, width_: Quantity,
    height_: Quantity) -> Quantity:
    if first_index_ < 0 or second_index_ < 0:
        raise ValueError("Indexes should not be negative")
    result_velocity_expr = solve(law, critical_wavelength, dict=True)[0][critical_wavelength]
    result_expr = result_velocity_expr.subs({
        first_index: first_index_,
        second_index: second_index_,
        width: width_,
        height: height_
    })
    return Quantity(result_expr)
