from sympy import (Eq, solve, log)
from symplyphysics import (
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
## Thermal effect of reaction is change of enthalpy of the system.
## The standard state is the state at a temperature of 298 Kelvin and a total pressure of 1 atmosphere,
## as well as at a fixed composition of the system.

## Law is: G = H - T * S - Cp * T * (ln(T / 298) + (298 / T) - 1), where
## G - change of isobaric potential of reaction,
## H - thermal effect of reaction,
## T - temperature,
## S - change of standart_change_entropy,
## Cp - change of heat capacity.

# Conditions:
## - we neglect the temperature dependence of the heat capacities;
## - the process is isobaric-isothermal.

standart_change_isobaric_potential = Symbol("standart_change_isobaric_potential", units.energy / units.amount_of_substance)

standart_thermal_effect = Symbol("standart_thermal_effect", units.energy / units.amount_of_substance)
standart_change_entropy = Symbol("standart_change_entropy", units.energy / units.amount_of_substance / units.temperature)
temperature = Symbol("temperature", units.temperature)
standart_change_heat_capacity = Symbol("standart_change_heat_capacity", units.energy / units.amount_of_substance / units.temperature)
standard_temperature = Quantity(298 * units.kelvin)

law = Eq(standart_change_isobaric_potential,
         standart_thermal_effect - temperature * standart_change_entropy - standart_change_heat_capacity * temperature * (log(temperature / standard_temperature) + (standard_temperature / temperature) - 1))


def print_law() -> str:
    return print_expression(law)


@validate_input(standart_thermal_effect_=standart_thermal_effect,
    standart_change_entropy_=standart_change_entropy,
    temperature_=temperature,
    standart_change_heat_capacity_=standart_change_heat_capacity)
@validate_output(standart_change_isobaric_potential)
def calculate_standart_change_isobaric_potential(standart_thermal_effect_: Quantity, standart_change_entropy_: Quantity,
    temperature_: Quantity, standart_change_heat_capacity_: Quantity) -> Quantity:
    result_expr = solve(law, standart_change_isobaric_potential, dict=True)[0][standart_change_isobaric_potential]
    result_expr = result_expr.subs({
        standart_thermal_effect: standart_thermal_effect_,
        standart_change_entropy: standart_change_entropy_,
        temperature: temperature_,
        standart_change_heat_capacity: standart_change_heat_capacity_,
    })
    return Quantity(result_expr)
