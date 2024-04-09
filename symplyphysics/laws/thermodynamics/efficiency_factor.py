from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
)
from symplyphysics.core.convert import convert_to_dimensionless

# Description
## The efficiency of a heat engine is the ratio of the useful energy used to the total amount of energy received by the system: eta = (Q_h - Q_r) / Q_h.
## Where:
## Q_h - the amount of heat transferred by the heater to the heat engine
## Q_r - the amount of heat transferred by the heat engine to the refrigerator
## eta - efficiency of the heat engine

heat_from_heater = Symbol("heat_from_heater", units.energy)
heat_to_refrigerator = Symbol("heat_to_refrigerator", units.energy)
efficiency_factor = Symbol("efficiency_factor", dimensionless)

law = Eq(efficiency_factor, (heat_from_heater - heat_to_refrigerator) / heat_from_heater)


def print_law() -> str:
    return print_expression(law)


@validate_input(heat_from_heater_=heat_from_heater, heat_to_refrigerator_=heat_to_refrigerator)
@validate_output(efficiency_factor)
def calculate_efficiency_factor(heat_from_heater_: Quantity,
    heat_to_refrigerator_: Quantity) -> float:
    result_efficiency_factor = solve(law, efficiency_factor, dict=True)[0][efficiency_factor]
    result_expr = result_efficiency_factor.subs({
        heat_from_heater: heat_from_heater_,
        heat_to_refrigerator: heat_to_refrigerator_
    })
    result = Quantity(result_expr)
    return convert_to_dimensionless(Quantity(result))
