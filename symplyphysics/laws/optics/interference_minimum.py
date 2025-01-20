"""
Interference minimum
====================

If the travel difference is equal to an odd number of half-wavelengths, then at this point of the screen
there will be a minimum of intensity during interference.

**Conditions:**

#. The waves must be coherent.

**Links:**

#. `Physics LibreTexts <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_III_-_Optics_and_Modern_Physics_(OpenStax)/03%3A_Interference/3.02%3A_Young%27s_Double-Slit_Interference>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (Quantity, validate_input, validate_output, symbols)

travel_difference = symbols.optical_distance
"""
Optical travel difference between two coherent waves. See :symbols:`optical_distance`.
"""

number_minimum = symbols.whole_number
"""
A :symbols:`whole_number` of interference minima.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength` of the two waves.
"""

law = Eq(travel_difference, (2 * number_minimum + 1) * wavelength / 2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(wave_length_=wavelength, number_minimum_=number_minimum)
@validate_output(travel_difference)
def calculate_travel_difference(wave_length_: Quantity, number_minimum_: int) -> Quantity:
    solved = solve(law, travel_difference, dict=True)[0][travel_difference]
    result_expr = solved.subs({wavelength: wave_length_, number_minimum: number_minimum_})
    return Quantity(result_expr)
