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
## Thermal effect of reaction is change of enthalpy of the system.
## The standard state is the state at a temperature of 298 Kelvin and a total pressure of 1 atmosphere,
## as well as at a fixed composition of the system.

## Law is: G = H - T * S - Cp * T * (ln(T / 298) + (298 / T) - 1), where
## G - standard change of isobaric potential of reaction,
## H - standard thermal effect of reaction,
## T - temperature,
## S - standard change of entropy,
## Cp - standard change of heat capacity.

# Conditions:
## - we neglect the temperature dependence of the heat capacities;
## - the process is isobaric-isothermal.

# TODO: find link

standard_change_isobaric_potential = Symbol("standard_change_isobaric_potential",
    units.energy / units.amount_of_substance)

standard_thermal_effect = Symbol("standard_thermal_effect",
    units.energy / units.amount_of_substance)
standard_change_entropy = Symbol("standard_change_entropy",
    units.energy / units.amount_of_substance / units.temperature)
temperature = symbols.temperature
standard_change_heat_capacity = Symbol("standard_change_heat_capacity",
    units.energy / units.amount_of_substance / units.temperature)

law = Eq(
    standard_change_isobaric_potential, standard_thermal_effect -
    temperature * standard_change_entropy - standard_change_heat_capacity * temperature *
    (log(temperature / quantities.standard_laboratory_temperature) +
    (quantities.standard_laboratory_temperature / temperature) - 1))


def print_law() -> str:
    return print_expression(law)


@validate_input(standard_thermal_effect_=standard_thermal_effect,
    standard_change_entropy_=standard_change_entropy,
    temperature_=temperature,
    standard_change_heat_capacity_=standard_change_heat_capacity)
@validate_output(standard_change_isobaric_potential)
def calculate_standard_change_isobaric_potential(
        standard_thermal_effect_: Quantity, standard_change_entropy_: Quantity,
        temperature_: Quantity, standard_change_heat_capacity_: Quantity) -> Quantity:
    result_expr = solve(law, standard_change_isobaric_potential,
        dict=True)[0][standard_change_isobaric_potential]
    result_expr = result_expr.subs({
        standard_thermal_effect: standard_thermal_effect_,
        standard_change_entropy: standard_change_entropy_,
        temperature: temperature_,
        standard_change_heat_capacity: standard_change_heat_capacity_,
    })
    return Quantity(result_expr)
