from sympy import (Eq, solve)
from symplyphysics import (Symbol, print_expression, validate_input,
    validate_output, dimensionless, convert_to_float)

# Description
## There is a coupling parameter to describe the resonator and the load. The parameter is equal to the ratio of the
## resonator's own quality factor to the quality factor of the external circuit.

## Law is: g = Q0 / Qe, where
## g - coupling parameter,
## Q0 - resonator's own quality factor,
## Qe - quality factor of the external circuit.

coupling_parameter = Symbol("coupling_parameter", dimensionless)

resonator_quality_factor = Symbol("resonator_quality_factor", dimensionless)
external_circuit_quality_factor = Symbol("external_circuit_quality_factor", dimensionless)

law = Eq(coupling_parameter, resonator_quality_factor / external_circuit_quality_factor)


def print_law() -> str:
    return print_expression(law)


@validate_input(resonator_quality_factor_=resonator_quality_factor, external_circuit_quality_factor_=external_circuit_quality_factor)
@validate_output(coupling_parameter)
def calculate_coupling_parameter(resonator_quality_factor_: float,
    external_circuit_quality_factor_: float) -> float:
    result_expr = solve(law, coupling_parameter, dict=True)[0][coupling_parameter]
    result_expr = result_expr.subs({
        resonator_quality_factor: resonator_quality_factor_,
        external_circuit_quality_factor: external_circuit_quality_factor_,
    })
    return convert_to_float(result_expr)
