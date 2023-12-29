from sympy import Eq, solve
from sympy.physics.units import joule
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    LessEq
)


# Description
## This inequality states that the efficiency of any heat engine must be less than or equal to the maximum value of 1 - Q1/Q2
## Efficiency according to the Kelvin-Planck law: eta <= 1 - Q2/Q1
## Where:
## eta - machine efficiency
## Q1 - quantity of heat taken from the hot reservoir in joules (initial temperature)
## Q2 - quantity of heat given to the cold reservoir in joules (final temperature)

eta = Symbol("efficiency", 1 - (units.joule / units.joule))
temperature_start = Symbol("temperature_start", units.joule)
temperature_end = Symbol("temperature_end", units.joule)

law = LessEq(eta, 1 - (temperature_end / temperature_start))

def print_law() -> str:
    return print_expression(law)


@validate_input(temperature_start_=temperature_start, temperature_end_=temperature_end)
@validate_output(eta)
def calculate_eta(
    temperature_start_: Quantity,
    temperature_end_: Quantity,
) -> Quantity:
    solved = solve(law, eta, dict=True)[0][eta]
    result_expr = solved.subs(
        {temperature_start: temperature_start_, temperature_end: temperature_end_}
    )
    return Quantity(result_expr)
