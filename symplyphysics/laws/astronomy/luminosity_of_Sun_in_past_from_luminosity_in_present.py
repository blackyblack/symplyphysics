from sympy import (Eq, solve, S,)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    convert_to,
    dimensionless,
)

# Description
## The luminosity of the Sun in the past is related to the luminosity of the Sun in the present. Luminosity is
## indicated in units of solar luminosity. One unit is equal to the luminosity of the Sun at a given time.

## Law is: Lpast = L / (1 + 0.4 * (1 - (t / 4.6))), where
## Lpast - luminosity of the sun in past,
## L - luminosity of the sun in present,
## t - time in billions of years.

luminosity_past = Symbol("luminosity_past", dimensionless)

luminosity_present = Symbol("luminosity_present", dimensionless)
time = Symbol("time", dimensionless)

law = Eq(luminosity_past, luminosity_present / (1 + 0.4 * (1 - (time / 4.6))))


def print_law() -> str:
    return print_expression(law)


@validate_input(luminosity_present_=luminosity_present, time_=time)
@validate_output(luminosity_past)
def calculate_luminosity_past(luminosity_present_: float, time_: float) -> float:
    result_expr = solve(law, luminosity_past, dict=True)[0][luminosity_past]
    result_expr = result_expr.subs({
        luminosity_present: luminosity_present_,
        time: time_,
    })
    return float(convert_to(Quantity(result_expr), S.One).evalf())
