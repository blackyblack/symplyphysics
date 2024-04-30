from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    convert_to_float,
    dimensionless,
)

# Description
## The reflection coefficient is equal to the ratio of the reflected wave power to the incident wave power.
## Knowing the average power delivered to the load and the incident power, it is possible to calculate the absolute value reflection coefficient.
## It is also stated in the description of the law that G is the modulus of the reflection coefficient. But in fact, G can also be a reflection
## coefficient (a complex number). But since we want to calculate the value of G uniquely from the law, we will assume that G is the absolute value
## of the reflection coefficient.

## Law is: Pa / Pi = 1 - G^2, where
## Pa - the average power delivered to the load,
## Pi - the incident power,
## G - the absolute value reflection coefficient.

absolute_reflection_coefficient = Symbol("absolute_reflection_coefficient", dimensionless, real=True)
incident_power = Symbol("incident_power", units.power, real=True)
average_power = Symbol("average_power", units.power, real=True)

law = Eq(average_power / incident_power, 1 - absolute_reflection_coefficient**2)


@validate_input(incident_power_=incident_power, average_power_=average_power)
@validate_output(absolute_reflection_coefficient)
def calculate_absolute_reflection_coefficient(incident_power_: Quantity, average_power_: Quantity) -> float:
    if incident_power_.scale_factor < average_power_.scale_factor:
        raise ValueError("The incident_power must be greater than the average power")
    result_expr = solve(law, absolute_reflection_coefficient, dict=True)[0][absolute_reflection_coefficient]
    result_expr = result_expr.subs({
        incident_power: incident_power_,
        average_power: average_power_,
    })
    return convert_to_float(result_expr)
