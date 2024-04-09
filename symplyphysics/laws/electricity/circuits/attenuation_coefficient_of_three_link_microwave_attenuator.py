from sympy import Eq, solve, exp, acosh
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
)

# Description
## Microwave attenuators are used to attenuate the microwave signal. For a three-link T-type attenuator or a π-type
## attenuator, the signal attenuation coefficient is calculated from the resistor resistance ratio.

## Law is: N = exp(arcosh(1 + R1 / R2)), where
## N - attenuation coefficient of attenuator,
## R1 - resistance of the first resistor,
## R2 - resistance of the second resistor,
## arcosh - this is area hyperbolic cosine (https://en.wikipedia.org/wiki/Hyperbolic_functions).

attenuation_coefficient = Symbol("attenuation_coefficient", dimensionless)

first_resistance = Symbol("first_resistance", units.impedance)
second_resistance = Symbol("second_resistance", units.impedance)

law = Eq(attenuation_coefficient, exp(acosh(1 + first_resistance / second_resistance)))


def print_law() -> str:
    return print_expression(law)


@validate_input(first_resistance_=first_resistance, second_resistance_=second_resistance)
@validate_output(attenuation_coefficient)
def calculate_attenuation_coefficient(first_resistance_: Quantity,
    second_resistance_: Quantity) -> float:
    result_expr = solve(law, attenuation_coefficient, dict=True)[0][attenuation_coefficient]
    result_expr = result_expr.subs({
        first_resistance: first_resistance_,
        second_resistance: second_resistance_,
    })
    return convert_to_float(Quantity(result_expr))
