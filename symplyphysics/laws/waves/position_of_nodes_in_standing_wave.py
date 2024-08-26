"""
Position of nodes in standing wave
==================================

In a standing wave, the locations of zero amplitude are called nodes. These locations, however,
are not arbitrary and are integer multiples of half the wavelength of the standing wave.
"""

from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import angular_wavenumber_is_inverse_wavelength as wavenumber_def
from symplyphysics.laws.waves import displacement_in_standing_wave as standing_wave_law

node_position = Symbol("node_position", units.length, real=True)
r"""
Position of :math:`m`-th node.

Symbol:
    :code:`x`
"""

integer_factor = Symbol("integer_factor", dimensionless, integer=True)
"""
An integer.

Symbol:
    :code:`m`
"""

wavelength = Symbol("wavelength", units.length, positive=True)
r"""
Wavelength of the standing wave.

Symbol:
    :code:`lambda`

Latex:
    :math:`\lambda`
"""

law = Eq(node_position, integer_factor * wavelength / 2)
r"""
:code:`x = (m / 2) * lambda`

Latex:
    .. math::
        x = \frac{m}{2} \lambda
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
