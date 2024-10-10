r"""
Wavelength of standing wave in string with fixed ends
=====================================================

If a standing wave occurs in a string with fixed ends, the value of the displacement
function :math:`q(x)` must be zero on both ends of the string. To satisfy this boundary
condition, only an integer number of half-waves can fit into the total length of the
string.

**Conditions:**

#. Boundary condition: :math:`q(0) = q(L) = 0`.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import (
    angular_wavenumber_is_inverse_wavelength as wavenumber_def,
)
from symplyphysics.laws.waves import displacement_in_standing_wave as standing_wave_law

integer_factor = symbols.positive_number
r"""
Positive integer, also called *harmonic number* of the :math:`n`-th harmonic.
See :symbols:`positive_number`.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength` of standing wave.
"""

string_length = symbols.length
"""
:symbols:`length` of the string.
"""

law = Eq(integer_factor * wavelength / 2, string_length)
"""
:laws:symbol::

:laws:latex::
"""

# Derive from boundary condition `u(L, t) = 0`

_wavenumber = wavenumber_def.definition.rhs.subs(
    wavenumber_def.wavelength,
    wavelength,
)

_standing_wave = standing_wave_law.law.rhs.subs(
    standing_wave_law.angular_wavenumber,
    _wavenumber,
)

_standing_wave_at_string_end = _standing_wave.subs(
    standing_wave_law.position,
    string_length,
)

_wavelength_solution = solve(_standing_wave_at_string_end, wavelength)[0]

# We can show that any wavelength that is a whole number of times shorter than the
# one we just got as the solution of the equation is also a resonant wavelength.

_wavelength_solution_set = _wavelength_solution / integer_factor

assert expr_equals(_standing_wave_at_string_end.subs(wavelength, _wavelength_solution_set), 0)

_wavelength_from_law = solve(law, wavelength)[0]

assert expr_equals(_wavelength_solution_set, _wavelength_from_law)


@validate_input(
    integer_factor_=integer_factor,
    string_length_=string_length,
)
@validate_output(wavelength)
def calculate_wavelength(
    integer_factor_: int,
    string_length_: Quantity,
) -> Quantity:
    result_expr = solve(law, wavelength)[0]
    result = result_expr.subs({
        integer_factor: integer_factor_,
        string_length: string_length_,
    })
    return Quantity(result)
