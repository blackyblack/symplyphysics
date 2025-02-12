"""
Standing wave ratio from reflection coefficient
===============================================

The standing wave ratio can be calculated by knowing the reflection coefficient of the
wave.

**Links:**

#. `Engineering LibreTexts <https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Electro-Optics/Book%3A_Electromagnetics_I_(Ellingson)/03%3A_Transmission_Lines/3.14%3A_Standing_Wave_Ratio>`__.

..
    TODO: fix file name
"""

from sympy import Eq, solve
from symplyphysics import validate_input, validate_output, convert_to_float, symbols

standing_wave_ratio = symbols.standing_wave_ratio
"""
:symbols:`standing_wave_ratio`.
"""

reflection_coefficient = symbols.reflection_coefficient
"""
:symbols:`reflection_coefficient`.
"""

law = Eq(standing_wave_ratio,
    (1 + abs(reflection_coefficient)) / (1 - abs(reflection_coefficient)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(reflection_coefficient_module_=reflection_coefficient)
@validate_output(standing_wave_ratio)
def calculate_coefficient_standing_wave(reflection_coefficient_module_: float) -> float:
    if reflection_coefficient_module_ < 0:
        raise ValueError("The modulus of the reflection coefficient cannot be less than zero")
    result_expr = solve(law, standing_wave_ratio, dict=True)[0][standing_wave_ratio]
    result_expr = result_expr.subs({
        reflection_coefficient: reflection_coefficient_module_,
    })
    return convert_to_float(result_expr)
