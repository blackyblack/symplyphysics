from sympy import Eq, solve
from symplyphysics import (Symbol, print_expression, validate_input, validate_output, dimensionless,
    convert_to_float)

# Description
## Standing wave coefficient is the ratio of the highest voltage in the transmission line to the lowest voltage. It is a measure of matching the load
## with the transmission line. The coefficient in the transmission line does not depend on the internal resistance of the electromagnetic
## wave source and (in the case of a linear load) on the power of the generator.
## Also, the standing wave coefficient can be calculated by knowing the reflection coefficient module.

## Law is: K = (1 + G) / (1 - G), where
## K - standing wave coefficient,
## G - reflection coefficient module.

coefficient_standing_wave = Symbol("coefficient_standing_wave", dimensionless)

reflection_coefficient_module = Symbol("reflection_coefficient_module", dimensionless)

law = Eq(coefficient_standing_wave,
    (1 + reflection_coefficient_module) / (1 - reflection_coefficient_module))


def print_law() -> str:
    return print_expression(law)


@validate_input(reflection_coefficient_module_=reflection_coefficient_module)
@validate_output(coefficient_standing_wave)
def calculate_coefficient_standing_wave(reflection_coefficient_module_: float) -> float:
    if reflection_coefficient_module_ < 0:
        raise ValueError("The modulus of the reflection coefficient cannot be less than zero")
    result_expr = solve(law, coefficient_standing_wave, dict=True)[0][coefficient_standing_wave]
    result_expr = result_expr.subs({
        reflection_coefficient_module: reflection_coefficient_module_,
    })
    return convert_to_float(result_expr)
