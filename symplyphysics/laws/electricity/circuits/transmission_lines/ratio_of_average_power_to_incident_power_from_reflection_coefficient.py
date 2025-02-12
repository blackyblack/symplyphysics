"""
Reflection coefficient from ratio of average power to incident power
====================================================================

Knowing the average power delivered to the load and the incident power, it is possible
to calculate the absolute value of the reflection coefficient.

**Links:**

#. `Engineering LibreTexts, formula 3.20.4 <https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Electro-Optics/Book%3A_Electromagnetics_I_(Ellingson)/03%3A_Transmission_Lines/3.20%3A_Power_Flow_on_Transmission_Lines>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

reflection_coefficient = symbols.reflection_coefficient
"""
Complex-valued :symbols:`reflection_coefficient`.
"""

incident_power = clone_as_symbol(symbols.power, display_symbol="P_incident", display_latex="P_\\text{incident}")
"""
Incident :symbols:`power`.
"""

average_power = clone_as_symbol(symbols.power, display_symbol="avg(P)", display_latex="\\langle P \\rangle")
"""
Average :symbols:`power` delivered to the load.
"""

law = Eq(average_power / incident_power, 1 - abs(reflection_coefficient)**2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(incident_power_=incident_power, average_power_=average_power)
@validate_output(reflection_coefficient)
def calculate_absolute_reflection_coefficient(incident_power_: Quantity,
    average_power_: Quantity) -> float:
    if incident_power_.scale_factor < average_power_.scale_factor:
        raise ValueError("The incident_power must be greater than the average power")
    abs_gamma = clone_as_symbol(reflection_coefficient, positive=True)
    subs_law = law.subs(abs(reflection_coefficient), abs_gamma)
    result_expr = solve(subs_law, abs_gamma)[1]
    result_expr = result_expr.subs({
        incident_power: incident_power_,
        average_power: average_power_,
    })
    return convert_to_float(result_expr)
