"""
Phase shift between two points
==============================

The phase shift of wave oscillations between two points is proportional to the ratio
of the distance between the points and the wavelength of the wave.

**Links:**

#. `BYJY's <https://byjus.com/physics/relation-between-phase-difference-and-path-difference/>`__.
"""

from sympy import Eq, solve, pi
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
)

phase_shift = symbols.phase_shift
"""
:symbols:`phase_shift` between the two points.
"""

distance = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` between the two points.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength` of the wave.
"""

law = Eq(phase_shift, 2 * pi * distance / wavelength)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(distance_between_points_=distance, wavelength_=wavelength)
@validate_output(phase_shift)
def calculate_phase_difference(distance_between_points_: Quantity, wavelength_: Quantity) -> float:
    result_expr = solve(law, phase_shift, dict=True)[0][phase_shift]
    result_expr = result_expr.subs({
        distance: distance_between_points_,
        wavelength: wavelength_
    }).doit()
    return convert_to_float(result_expr)
