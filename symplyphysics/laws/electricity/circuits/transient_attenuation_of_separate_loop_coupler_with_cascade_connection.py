from sympy import (
    Eq,
    solve,
    sin,
    asin,
    log
)
from symplyphysics import (Symbol, print_expression, validate_input,
    validate_output, dimensionless, convert_to_float)

# Description
## The multistage connection of the couplers allows you to expand the working band and realize less transient attenuation.
## Transient attenuation is a part of the power passing into the auxiliary channel of the coupler.

## Law is: Sn = 20 * lg(sin(n * asin(10^(S / 20)))), where
## Sn - transient attenuation of one coupler from the cascade (measured in decibels),
## n - the number of couplers in the cascade,
## S - transient attenuation of the cascade (measured in decibels).

attenuation_of_coupler = Symbol("attenuation_of_coupler", dimensionless)

attenuation_of_cascade = Symbol("attenuation_of_cascade", dimensionless)
number_of_couplers = Symbol("number_of_couplers", dimensionless)

law = Eq(attenuation_of_coupler, 20 * log(sin(number_of_couplers * asin(10**(attenuation_of_cascade / 20))), 10))


def print_law() -> str:
    return print_expression(law)


@validate_input(attenuation_of_cascade_=attenuation_of_cascade,
    number_of_couplers_=number_of_couplers)
@validate_output(attenuation_of_coupler)
def calculate_attenuation_of_coupler(attenuation_of_cascade_: float,
    number_of_couplers_: float) -> float:
    result_expr = solve(law, attenuation_of_coupler, dict=True)[0][attenuation_of_coupler]
    result_expr = result_expr.subs({
        attenuation_of_cascade: attenuation_of_cascade_,
        number_of_couplers: number_of_couplers_,
    })
    return convert_to_float(result_expr)
