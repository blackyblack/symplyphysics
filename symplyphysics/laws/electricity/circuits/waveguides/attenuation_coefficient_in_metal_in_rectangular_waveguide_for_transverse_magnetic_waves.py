"""
Attenuation coefficient in metal in rectangular waveguide for transverse magnetic waves
=======================================================================================

A coaxial waveguide is an electrical cable consisting of a central conductor and a
shield arranged coaxially and separated by an insulating material or an air gap. It is
used to transmit radio frequency electrical signals. The specific resistance of a
coaxial waveguide depends on the radius of the outer conductor and the radius of the
inner conductor, as well as on the relative permeability of the insulator material, 
frequency of signal and specific conductivity of conductor.

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
    SymbolNew,
)

attenuation_coefficient = symbols.attenuation_coefficient
"""
:symbols:`attenuation_coefficient` in metal.
"""

surface_resistance = clone_as_symbol(symbols.electrical_resistance, subscript="\\text{s}")
"""
:symbols:`electrical_resistance` of the surface.
"""

first_index = SymbolNew("m", dimensionless)
"""
The first index shows how many half-wavelengths fit across the width of the cross
section.
"""

second_index = SymbolNew("n", dimensionless)
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

medium_resistance = symbols.electrical_resistance
"""
:symbols:`electrical_resistance` of the medium filling the waveguide.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength` of the signal.
"""

critical_wavelength = clone_as_symbol(symbols.wavelength, subscript="\\text{c}")
"""
Critical :symbols:`wavelength` of the system. See :ref:`Critical wavelength of waveguide`.
"""

_reduced_resistance = surface_resistance / medium_resistance
_reduced_dimension = width / height
_reduced_wavelength = wavelength / (2 * critical_wavelength)

law = Eq(
    attenuation_coefficient,
    (2 * _reduced_resistance * (second_index**2 * _reduced_dimension**3 + first_index**2)) 
    / (
        width
        * sqrt(1 - _reduced_wavelength**2)
        * (second_index**2 * _reduced_dimension**2 + first_index**2)
    ),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(surface_resistance_=surface_resistance,
    first_index_=first_index,
    second_index_=second_index,
    width_=width,
    height_=height,
    resistance_of_medium_=medium_resistance,
    signal_wavelength_=wavelength,
    critical_wavelength_=critical_wavelength)
@validate_output(attenuation_coefficient)
def calculate_attenuation_coefficient(surface_resistance_: Quantity, first_index_: float,
    second_index_: float, width_: Quantity, height_: Quantity, resistance_of_medium_: Quantity,
    signal_wavelength_: Quantity, critical_wavelength_: Quantity) -> Quantity:
    # pylint: disable=too-many-arguments, too-many-positional-arguments
    result_expr = solve(law, attenuation_coefficient, dict=True)[0][attenuation_coefficient]
    result_expr = result_expr.subs({
        surface_resistance: surface_resistance_,
        first_index: first_index_,
        second_index: second_index_,
        width: width_,
        height: height_,
        medium_resistance: resistance_of_medium_,
        wavelength: signal_wavelength_,
        critical_wavelength: critical_wavelength_,
    })
    return Quantity(result_expr)
