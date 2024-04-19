from sympy import (Eq, solve)
from symplyphysics import (Symbol, print_expression, validate_input,
    validate_output, dimensionless, convert_to_float)

# Description
## The loaded quality factor of the resonator depends on the own quality factor of the resonator and
## the quality factor of the external circuit connected to the resonator.

## Law is: 1 / Ql = 1 / Q0 + 1 / Qe, where
## Ql - loaded resonator quality factor,
## Q0 - resonator's own quality factor,
## Qe - quality factor of the external circuit.

loaded_resonator_quality_factor = Symbol("loaded_resonator_quality_factor", dimensionless)

resonator_quality_factor = Symbol("resonator_quality_factor", dimensionless)
external_circuit_quality_factor = Symbol("external_circuit_quality_factor", dimensionless)

law = Eq(1 / loaded_resonator_quality_factor, 1 / resonator_quality_factor + 1 / external_circuit_quality_factor)


def print_law() -> str:
    return print_expression(law)


@validate_input(resonator_quality_factor_=resonator_quality_factor, external_circuit_quality_factor_=external_circuit_quality_factor)
@validate_output(loaded_resonator_quality_factor)
def calculate_quality_factor(resonator_quality_factor_: float,
    external_circuit_quality_factor_: float) -> float:
    result_expr = solve(law, loaded_resonator_quality_factor, dict=True)[0][loaded_resonator_quality_factor]
    result_expr = result_expr.subs({
        resonator_quality_factor: resonator_quality_factor_,
        external_circuit_quality_factor: external_circuit_quality_factor_,
    })
    return convert_to_float(result_expr)
