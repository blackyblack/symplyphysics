"""
Standing wave ratio from ratio of average power to incident power
=================================================================

Knowing the average power delivered to the load and the incident power, it is possible
to calculate the standing wave ratio.

..
    TODO: find link
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

standing_wave_ratio = symbols.standing_wave_ratio
"""
:symbols:`standing_wave_ratio`.
"""

incident_power = clone_as_symbol(symbols.power, display_symbol="P_incident", display_latex="P_\\text{incident}")
"""
Incident :symbols:`power`.
"""

average_power = clone_as_symbol(symbols.power, display_symbol="avg(P)", display_latex="\\langle P \\rangle")
"""
Average :symbols:`power` delivered to the load.
"""

law = Eq(average_power / incident_power,
    4 * standing_wave_ratio / (standing_wave_ratio + 1)**2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(incident_power_=incident_power, average_power_=average_power)
@validate_output(standing_wave_ratio)
def calculate_standing_wave_coefficient(incident_power_: Quantity,
    average_power_: Quantity) -> float:
    if incident_power_.scale_factor < average_power_.scale_factor:
        raise ValueError("The incident_power must be greater than the average power")
    result_expr = solve(law, standing_wave_ratio, dict=True)[1][standing_wave_ratio]
    result_expr = result_expr.subs({
        incident_power: incident_power_,
        average_power: average_power_,
    })
    return convert_to_float(result_expr)
