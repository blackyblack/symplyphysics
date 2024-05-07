from sympy import Eq, solve, log
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

## Description
## The microstrip line is a dielectric substrate on which a metal strip is applied.
## The inductance of the upper metal strip without taking into account the lower electrode can be
## calculated by knowing the geometric parameters of the strip.

## Law is: L = 2e-4 * l * (ln(l / (W + t)) + 1.193 + 0.2234 * ln((W + t) / l)), where
## L - inductance of the microstrip line strip,
## W - strip width of the microstrip line,
## l - strip length,
## t - strip thickness.

inductance = Symbol("inductance", units.inductance)

strip_thickness = Symbol("strip_thickness", units.length)
strip_length = Symbol("strip_length", units.length)
strip_width = Symbol("strip_width", units.length)

expression_1 = strip_length / (strip_width + strip_thickness)
constant_inductance = Quantity(2e-4 * units.inductance / units.length)

law = Eq(inductance, constant_inductance * strip_length * (log(expression_1) + 1.193 + 0.2235 * (1 / expression_1)))


def print_law() -> str:
    return print_expression(law)


@validate_input(strip_thickness_=strip_thickness,
    strip_length_=strip_length,
    strip_width_=strip_width)
@validate_output(inductance)
def calculate_inductance(strip_thickness_: Quantity, strip_length_: Quantity,
    strip_width_: Quantity) -> Quantity:
    result_expr = solve(law, inductance, dict=True)[0][inductance]
    result_expr = result_expr.subs({
        strip_thickness: strip_thickness_,
        strip_length: strip_length_,
        strip_width: strip_width_
    })
    return Quantity(result_expr)
