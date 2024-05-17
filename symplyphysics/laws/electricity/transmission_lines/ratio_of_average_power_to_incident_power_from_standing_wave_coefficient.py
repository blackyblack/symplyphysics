from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
)

# Description
## Standing wave coefficient is the ratio of the highest voltage in the transmission line to the lowest voltage. It is a measure of matching the load
## with the transmission line. The coefficient in the transmission line does not depend on the internal resistance of the electromagnetic
## wave source and (in the case of a linear load) on the power of the generator.
## Knowing the average power delivered to the load and the incident power, it is possible to calculate the standing wave coefficient.

## Law is: Pa / Pi = 4 * K / (K + 1)^2, where
## Pa - the average power delivered to the load,
## Pi - the incident power,
## K - the standing wave coefficient.

standing_wave_coefficient = Symbol("standing_wave_coefficient", dimensionless)
incident_power = Symbol("incident_power", units.power)
average_power = Symbol("average_power", units.power)

law = Eq(average_power / incident_power,
    4 * standing_wave_coefficient / (standing_wave_coefficient + 1)**2)


@validate_input(incident_power_=incident_power, average_power_=average_power)
@validate_output(standing_wave_coefficient)
def calculate_standing_wave_coefficient(incident_power_: Quantity,
    average_power_: Quantity) -> float:
    if incident_power_.scale_factor < average_power_.scale_factor:
        raise ValueError("The incident_power must be greater than the average power")
    result_expr = solve(law, standing_wave_coefficient, dict=True)[1][standing_wave_coefficient]
    result_expr = result_expr.subs({
        incident_power: incident_power_,
        average_power: average_power_,
    })
    return convert_to_float(result_expr)
