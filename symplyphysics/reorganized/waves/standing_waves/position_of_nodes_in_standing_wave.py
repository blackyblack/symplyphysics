"""
Position of nodes in standing wave
==================================

In a standing wave, the locations of zero amplitude are called nodes. These locations, however,
are not arbitrary and are integer multiples of half the wavelength of the standing wave.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Standing_wave#Standing_wave_on_an_infinite_length_string>`__.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import angular_wavenumber_is_inverse_wavelength as wavenumber_def
from symplyphysics.laws.waves import displacement_in_standing_wave as standing_wave_law

node_position = symbols.position
r"""
:symbols:`position` of :math:`m`-th node.
"""

integer_factor = symbols.whole_number
"""
An integer. See :symbols:`whole_number`.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength` of the standing wave.
"""

law = Eq(node_position, integer_factor * wavelength / 2)
"""
:laws:symbol::

:laws:latex::
"""

# Proving these are indeed locations of zero amplitude.
# Deriving it from the standing wave expression is impossible since `sympy` does not produce
# infinite solutions for trigonometric equations.

_standing_wave_expr = standing_wave_law.law.rhs

_angular_wavenumber = wavenumber_def.definition.rhs.subs(wavenumber_def.wavelength, wavelength)

_standing_wave_expr = _standing_wave_expr.subs(standing_wave_law.angular_wavenumber,
    _angular_wavenumber)

_node_amplitude = _standing_wave_expr.subs(standing_wave_law.position, law.rhs)

assert expr_equals(_node_amplitude, 0)


@validate_input(
    integer_factor_=integer_factor,
    wavelength_=wavelength,
)
@validate_output(node_position)
def calculate_node_position(
    integer_factor_: int,
    wavelength_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        integer_factor: integer_factor_,
        wavelength: wavelength_,
    })
    return Quantity(result)
