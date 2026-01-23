"""
Interference maximum
====================

If the difference in the optical path of two waves is equal to an integer number of wavelengths
(i.e. an even number of half-wavelengths) then an interference maximum is formed at the point
of superposition of these waves.

**Conditions:**

#. The waves must be coherent.

**Links:**

#. `Vogueindustry <https://vogueindustry.com/17289645-interference-patterns-maximum-and-minimum-conditions#menu-9>`__.
#. `Livelaptopspec <https://www.livelaptopspec.com/what-is-maximum-constructive-interference/>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (Quantity, validate_input, validate_output, symbols)

travel_difference = symbols.optical_distance
"""
Optical travel difference between two coherent waves. See :symbols:`optical_distance`.
"""

integer_factor = symbols.whole_number
"""
A :symbols:`whole_number` of interference maxima.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength` of each of the both waves.
"""

law = Eq(travel_difference, integer_factor * wavelength)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(wave_length_=wavelength, number_maximum_=integer_factor)
@validate_output(travel_difference)
def calculate_travel_difference(wave_length_: Quantity, number_maximum_: int) -> Quantity:
    solved = solve(law, travel_difference, dict=True)[0][travel_difference]
    result_expr = solved.subs({wavelength: wave_length_, integer_factor: number_maximum_})
    return Quantity(result_expr)
