from sympy import (Eq, solve)
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
## Thermal effect of reaction is enthalpy of the system.

## Law is: G = H - T * S, where
## G - isobaric potential of reaction,
## H - thermal effect of reaction,
## T - temperature,
## S - entropy.

# Conditions:
## - we neglect the temperature dependence of thermodynamic quantities;
## - the process is isobaric-isothermal.

isobaric_potential = Symbol("isobaric_potential", units.energy / units.amount_of_substance)

thermal_effect = Symbol("thermal_effect", units.energy / units.amount_of_substance)
entropy = Symbol("entropy", units.energy / units.amount_of_substance / units.temperature)
temperature = Symbol("temperature", units.temperature)

law = Eq(isobaric_potential, thermal_effect - temperature * entropy)


def print_law() -> str:
    return print_expression(law)


@validate_input(thermal_effect_=thermal_effect,
    entropy_=entropy,
    temperature_=temperature)
@validate_output(isobaric_potential)
def calculate_isobaric_potential(thermal_effect_: Quantity, entropy_: Quantity,
    temperature_: Quantity) -> Quantity:
    result_expr = solve(law, isobaric_potential, dict=True)[0][isobaric_potential]
    result_expr = result_expr.subs({
        thermal_effect: thermal_effect_,
        entropy: entropy_,
        temperature: temperature_
    })
    return Quantity(result_expr)
