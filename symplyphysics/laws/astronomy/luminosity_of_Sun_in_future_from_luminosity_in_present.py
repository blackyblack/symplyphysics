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
## The luminosity of the Sun in the future is related to the luminosity of the Sun in the present. Luminosity is
## indicated in units of solar luminosity. One unit is equal to the luminosity of the Sun at a given time.

## Law is: Lf = L * ((5.59 / t) - 1.41 + 0.26 * t), where
## Lf - luminosity of the sun in future,
## L - luminosity of the sun in present,
## t - time in billions of years.

luminosity_future = Symbol("luminosity_future", dimensionless)

luminosity_present = Symbol("luminosity_present", dimensionless)
time = Symbol("time", units.time)

time_constant = Quantity(1e9 * units.common_year)

law = Eq(luminosity_future, luminosity_present * ((5.59 / (time / time_constant)) - 1.41 + 0.26 * (time / time_constant)))


def print_law() -> str:
    return print_expression(law)


@validate_input(luminosity_present_=luminosity_present, time_=time)
@validate_output(luminosity_future)
def calculate_luminosity_future(luminosity_present_: float, time_: Quantity) -> float:
    result_expr = solve(law, luminosity_future, dict=True)[0][luminosity_future]
    result_expr = result_expr.subs({
        luminosity_present: luminosity_present_,
        time: time_,
    })
    return float(convert_to(Quantity(result_expr), S.One).evalf())
