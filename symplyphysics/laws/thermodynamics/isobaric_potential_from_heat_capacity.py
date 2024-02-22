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

## Law is: G = H - T * S - Cp * T * (ln(T / 298) + (298 / T) - 1), where
## G - change of isobaric potential of reaction,
## H - thermal effect of reaction,
## T - temperature,
## S - change of entropy,
## Cp - change of heat capacity.

# Conditions:
## - we neglect the temperature dependence of the heat capacities;
## - the process is isobaric-isothermal.

isobaric_potential = Symbol("isobaric_potential", units.energy / units.amount_of_substance)

thermal_effect = Symbol("thermal_effect", units.energy / units.amount_of_substance)
entropy = Symbol("entropy", units.energy / units.amount_of_substance / units.temperature)
temperature = Symbol("temperature", units.temperature)
heat_capacity = Symbol("heat_capacity", units.energy / units.amount_of_substance / units.temperature)
standard_temperature = Quantity(298 * units.kelvin)

law = Eq(isobaric_potential,
         thermal_effect - temperature * entropy - heat_capacity * temperature * (log(temperature / standard_temperature) + (standard_temperature / temperature) - 1))


def print_law() -> str:
    return print_expression(law)


@validate_input(thermal_effect_=thermal_effect,
    entropy_=entropy,
    temperature_=temperature,
    heat_capacity_=heat_capacity)
@validate_output(isobaric_potential)
def calculate_isobaric_potential(thermal_effect_: Quantity, entropy_: Quantity,
    temperature_: Quantity, heat_capacity_: Quantity) -> Quantity:
    result_expr = solve(law, isobaric_potential, dict=True)[0][isobaric_potential]
    result_expr = result_expr.subs({
        thermal_effect: thermal_effect_,
        entropy: entropy_,
        temperature: temperature_,
        heat_capacity: heat_capacity_,
    })
    return Quantity(result_expr)
