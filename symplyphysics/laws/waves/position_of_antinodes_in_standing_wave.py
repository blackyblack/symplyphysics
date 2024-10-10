"""
Position of antinodes in standing wave
======================================

In a standing wave, the locations of maximum amplitude are called antinodes. These locations are
not arbitrary, however, and are integer multiples of half the wavelength shifted by a quarter
of the wavelength.

**Notes:**

#. *Amplitude* refers to the absolute value of displacement from rest, therefore positions of
   maximum amplitude are such positions where the displacement is maximum or minimum .
"""

from sympy import Eq, S
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import angular_wavenumber_is_inverse_wavelength as wavenumber_def
from symplyphysics.laws.waves import displacement_in_standing_wave as standing_wave_law

antinode_position = symbols.position
r"""
:symbols:`position` of :math:`m`-th antinode.
"""

integer_factor = symbols.whole_number
"""
An integer. See :symbols:`whole_number`.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength` of the standing wave.
"""

law = Eq(antinode_position, (integer_factor + S.One / 2) * wavelength / 2)
"""
:laws:symbol::

:laws:latex::
"""

# Proving these are indeed locations of maximum amplitude.
# Deriving it from the standing wave expression is impossible since `sympy` does not produce
# infinite solutions for trigonometric equations.

_standing_wave_expr = standing_wave_law.law.rhs

_angular_wavenumber = wavenumber_def.definition.rhs.subs(wavenumber_def.wavelength, wavelength)

_standing_wave_expr = _standing_wave_expr.subs(standing_wave_law.angular_wavenumber,
    _angular_wavenumber)
_standing_wave_spatial_derivative = _standing_wave_expr.diff(standing_wave_law.position)

_standing_wave_spatial_derivative_at_antinodes = _standing_wave_spatial_derivative.subs(
    standing_wave_law.position, law.rhs)

assert expr_equals(_standing_wave_spatial_derivative_at_antinodes, 0)

# Since the spatial derivative of the wave is zero in antinodes, according to calculus
# these are the points of extrema of the wave, hence they are locations of maximum and
# minimum displacement of the wave, hence they are locations of maximum amplitude.


@validate_input(
    integer_factor_=integer_factor,
    wavelength_=wavelength,
)
@validate_output(antinode_position)
def calculate_antinode_position(
    integer_factor_: int,
    wavelength_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        integer_factor: integer_factor_,
        wavelength: wavelength_,
    })
    return Quantity(result)
