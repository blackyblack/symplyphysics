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
## The luminosity of the Sun in the past is related to the luminosity of the Sun in the present. Luminosity is
## indicated in units of solar luminosity. One unit is equal to the luminosity of the Sun at a given time.
## There are two different formulas to calculate luminosity in the past and in the future. It is because it is
## non linear and requires different approximations for past and future.

## Law is: Lpast = L / (1 + 0.4 * (1 - (t / 4.6))), where
## Lpast - luminosity of the sun in past,
## L - luminosity of the sun in present,
## t - time in billions of years.

luminosity_past = Symbol("luminosity_past", dimensionless)

luminosity_present = Symbol("luminosity_present", dimensionless)
time = Symbol("time", units.time)

one_billion_years = Quantity(1e9 * units.common_year)

law = Eq(luminosity_past, luminosity_present / (1 + 0.4 * (1 - ((time / one_billion_years) / 4.6))))


def print_law() -> str:
    return print_expression(law)


@validate_input(luminosity_present_=luminosity_present, time_=time)
@validate_output(luminosity_past)
def calculate_luminosity_past(luminosity_present_: float, time_: Quantity) -> float:
    result_expr = solve(law, luminosity_past, dict=True)[0][luminosity_past]
    result_expr = result_expr.subs({
        luminosity_present: luminosity_present_,
        time: time_,
    })
    return convert_to_dimensionless(Quantity(result_expr))
