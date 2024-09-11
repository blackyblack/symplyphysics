from sympy import (Eq, solve, log)
from symplyphysics import (
    quantities,
    symbols,
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The isobaric potential of a reaction is a value whose change during a chemical reaction is equal to the change in the internal
## energy of the system. The isobaric potential shows how much of the total internal energy of the system can be used for chemical
## transformations.
## Thermal effect of reaction is enthalpy of the system.
## The heat capacity coefficients are tabular values for the reaction. They are used to express the dependence of heat capacity on temperature.

## Law is: G = H - T * S - T * (da * (ln(298 / T) + (298 / T) - 1) + db * ((T / 2) + (298^2 / (2 * T)) - 298) + dc * ((T^(-2) / 2) + (298^(-1)/-T) - (298^(-2) / -2))), where
## G - isobaric potential of reaction,
## H - thermal effect of reaction,
## T - temperature,
## S - entropy,
## da - the first tabular coefficient of heat capacity,
## db - the second tabular coefficient of heat capacity,
## dc - the third tabular coefficient of heat capacity.

# Conditions:
## - we take into account the dependence of heat capacity on temperature according to the Temkin-Schwarzman formula;
## - the process is isobaric-isothermal.

isobaric_potential = Symbol("isobaric_potential", units.energy / units.amount_of_substance)

thermal_effect = Symbol("thermal_effect", units.energy / units.amount_of_substance)
entropy = Symbol("entropy", units.energy / units.amount_of_substance / units.temperature)
temperature = symbols.temperature
heat_capacity = Symbol("heat_capacity",
    units.energy / units.amount_of_substance / units.temperature)
coefficient_capacity_1 = Symbol("coefficient_capacity_1",
    units.energy / units.amount_of_substance / units.temperature)
coefficient_capacity_2 = Symbol("coefficient_capacity_2",
    units.energy / units.amount_of_substance / units.temperature**2)
coefficient_capacity_3 = Symbol("coefficient_capacity_3",
    units.energy * units.temperature / units.amount_of_substance)

law = Eq(
    isobaric_potential, thermal_effect - temperature * entropy - temperature *
    (coefficient_capacity_1 * (log(quantities.standard_laboratory_temperature / temperature) +
    (quantities.standard_laboratory_temperature / temperature) - 1) + coefficient_capacity_2 *
    ((temperature / 2) + (quantities.standard_laboratory_temperature**2 /
    (2 * temperature)) - quantities.standard_laboratory_temperature) + coefficient_capacity_3 *
    ((temperature**(-2) / 2) + (quantities.standard_laboratory_temperature**(-1) / -temperature) -
    (quantities.standard_laboratory_temperature**(-2) / -2))))


def print_law() -> str:
    return print_expression(law)


#pylint: disable=too-many-arguments
@validate_input(
    thermal_effect_=thermal_effect,
    entropy_=entropy,
    temperature_=temperature,
    heat_capacity_=heat_capacity,
    coefficient_capacity_1_=coefficient_capacity_1,
    coefficient_capacity_2_=coefficient_capacity_2,
    coefficient_capacity_3_=coefficient_capacity_3,
)
@validate_output(isobaric_potential)
def calculate_isobaric_potential(thermal_effect_: Quantity, entropy_: Quantity,
    temperature_: Quantity, heat_capacity_: Quantity, coefficient_capacity_1_: Quantity,
    coefficient_capacity_2_: Quantity, coefficient_capacity_3_: Quantity) -> Quantity:
    result_expr = solve(law, isobaric_potential, dict=True)[0][isobaric_potential]
    result_expr = result_expr.subs({
        thermal_effect: thermal_effect_,
        entropy: entropy_,
        temperature: temperature_,
        heat_capacity: heat_capacity_,
        coefficient_capacity_1: coefficient_capacity_1_,
        coefficient_capacity_2: coefficient_capacity_2_,
        coefficient_capacity_3: coefficient_capacity_3_,
    })
    return Quantity(result_expr)
