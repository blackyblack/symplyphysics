"""
Maximum electric field strength of main wave in rectangular waveguide
=====================================================================

A rectangular waveguide is a rectangular metal waveguide capable of supporting waves
propagating along it. The main wave is a transverse electric wave with the first index
equal to :math:`1` and the second index equal to :math:`0`.

**Notes:**

#. *Horizontal* dimension refers to the dimension of the cross section of the waveguide
   that is orthogonal to its central axis. See image for reference.

.. image:: https://www.electronics-notes.com/images/waveguide-rectangular-te-modes-01.svg
    :width: 400px
    :align: center

**Notation:**

#. :quantity_notation:`vacuum_impedance`.

**Conditions:**

#. The wave propagating in the waveguide must be the main wave.
#. The waveguide must be rectangular.

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
    quantities,
)

maximum_electric_field_strength = symbols.electric_field_strength
"""
Maximum :symbols:`electric_field_strength` in the waveguide.
"""

relative_permittivity = symbols.relative_permittivity
"""
:symbols:`relative_permittivity` of the insulator.
"""

width = clone_as_symbol(symbols.length, display_symbol="a", display_latex="a")
"""
Horizontal dimension of the waveguide. See :symbols:`length`.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength` of the signal.
"""

magnetic_field_strength = symbols.magnetic_field_strength
"""
:symbols:`magnetic_field_strength`.
"""

law = Eq(
    maximum_electric_field_strength, 2 * quantities.vacuum_impedance * width * magnetic_field_strength /
    (wavelength * sqrt(relative_permittivity)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permittivity_=relative_permittivity,
    waveguide_width_=width,
    wavelength_=wavelength,
    magnetic_intensity_=magnetic_field_strength)
@validate_output(maximum_electric_field_strength)
def calculate_maximum_electric_intensity(relative_permittivity_: float, waveguide_width_: Quantity,
    wavelength_: Quantity, magnetic_intensity_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, maximum_electric_field_strength,
        dict=True)[0][maximum_electric_field_strength]
    result_expr = result_velocity_expr.subs({
        relative_permittivity: relative_permittivity_,
        width: waveguide_width_,
        wavelength: wavelength_,
        magnetic_field_strength: magnetic_intensity_
    })
    return Quantity(result_expr)
