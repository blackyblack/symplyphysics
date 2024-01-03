from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression,
                           validate_input, validate_output, dimensionless)

# Description
## Efficiency is a characteristic of the efficiency_law of a system (device, machine) in relation to the conversion or transfer of energy.
## It is determined by the ratio of the useful energy used to the total amount of energy received by the system: $\eta = (Q_h - Q_r) / Q_h$.
## Where:
## Q_h - the amount of heat transferred by the heater to the heat engine
## Q_r - the amount of heat transferred by the heat engine to the refrigerator
## \eta - efficiency_law of the heat engine

the_amount_of_heat_from_heater = Symbol("the_amount_of_heat_from_heater", units.joule)
the_amount_of_heat_to_refrigerator = Symbol("the_amount_of_heat_to_refrigerator", units.joule)
efficiency_factor = Symbol("efficiency_factor", dimensionless)

law = Eq(efficiency_factor, (the_amount_of_heat_from_heater - the_amount_of_heat_to_refrigerator) / the_amount_of_heat_from_heater)


def print_law() -> str:
    return print_expression(law)


@validate_input(the_amount_of_heat_from_heater_=the_amount_of_heat_from_heater,
                the_amount_of_heat_to_refrigerator_=the_amount_of_heat_to_refrigerator)
@validate_output(efficiency_factor)
def calculate_efficiency_factor(the_amount_of_heat_from_heater_: Quantity,
                                the_amount_of_heat_to_refrigerator_: Quantity) -> Quantity:
    result_efficiency_factor = solve(law, efficiency_factor, dict=True)[0][efficiency_factor]
    result_expr = result_efficiency_factor.subs({
        the_amount_of_heat_from_heater: the_amount_of_heat_from_heater_,
        the_amount_of_heat_to_refrigerator: the_amount_of_heat_to_refrigerator_
    })
    return Quantity(result_expr)
