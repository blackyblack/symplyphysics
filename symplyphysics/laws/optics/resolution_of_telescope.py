"""
Resolution of telescope
=======================

The resolution of the telescope is the minimum angular distance between point objects that can be
distinguished separately in the telescope.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Angular_resolution#The_Rayleigh_criterion>`__.
"""

from sympy import Eq, solve
from symplyphysics import units, Quantity, validate_input, validate_output, symbols

resolution = symbols.angular_resolution
"""
:symbols:`angular_resolution` of the telescope.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength` of light.
"""

lens_diameter = symbols.diameter
"""
:symbols:`diameter` of the lens.
"""

constant_rad = Quantity(1.22 * units.radian, display_symbol="A")
"""
Constant factor derived from a calculation of the position of the first dark circular ring surrounding the central Airy disc of the diffraction pattern. A :math:`\\approx` 1.22 radian.
"""

law = Eq(resolution, constant_rad * (wavelength / lens_diameter))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(wavelength_=wavelength, lens_diameter_=lens_diameter)
@validate_output(resolution)
def calculate_resolution(wavelength_: Quantity, lens_diameter_: Quantity) -> Quantity:
    result_expr = solve(law, resolution, dict=True)[0][resolution]
    result_expr = result_expr.subs({
        wavelength: wavelength_,
        lens_diameter: lens_diameter_,
    })
    return Quantity(result_expr)
