from sympy import (Eq, solve)
from symplyphysics import (units, Symbol, print_expression, validate_input,
    validate_output, dimensionless, convert_to_float)

# Description
## There is a coupling parameter to describe the resonator and the load. The parameter is equal to the ratio of the
## resonator's resistance to the load resistance.

## Law is: g = R0 / Rl, where
## g - coupling parameter,
## R0 - resonator's resistance,
## Rl - load resistance.

coupling_parameter = Symbol("coupling_parameter", dimensionless)

resonator_resistance = Symbol("resonator_resistance", units.impedance)
load_resistance = Symbol("load_resistance", units.impedance)

law = Eq(coupling_parameter, resonator_resistance / load_resistance)


def print_law() -> str:
    return print_expression(law)


@validate_input(resonator_resistance_=resonator_resistance, load_resistance_=load_resistance)
@validate_output(coupling_parameter)
def calculate_coupling_parameter(resonator_resistance_: float,
    load_resistance_: float) -> float:
    result_expr = solve(law, coupling_parameter, dict=True)[0][coupling_parameter]
    result_expr = result_expr.subs({
        resonator_resistance: resonator_resistance_,
        load_resistance: load_resistance_,
    })
    return convert_to_float(result_expr)
