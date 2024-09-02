"""
Phase shift between two points
==============================

The phase shift of wave oscillations between two points is proportional to the ratio
of the distance between the points and the wavelength of the wave.
"""

from sympy import Eq, solve, pi
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    angle_type,
    convert_to_float,
)

phase_shift = Symbol("phase_shift", angle_type)
r"""
Phase shift between the two points.

Symbol:
    :code:`phi`

Latex:
    :math:`\varphi`
"""

distance = Symbol("distance", units.length)
"""
Distance between the two points.

Symbol:
    :code:`d`
"""

wavelength = Symbol("wavelength", units.length)
r"""
Wavelength of the wave.

Symbol:
    :code:`lambda`

Latex:
    :math:`\lambda`
"""

law = Eq(phase_shift, 2 * pi * distance / wavelength)
r"""
:code:`phi = 2 * pi * d / lambda`

Latex:
    .. math::
        \varphi = 2 \pi \frac{d}{\lambda}
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
