from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## Kepler's third law is as follows. The squares of the periods of rotation of the planets around the Sun are referred to
## as cubes of the large semi-axes of the orbits of the planets.

## Law is: T1^2 / T2^2 = a1^3 / a2^3, where
## T1 - rotation period of the first planet,
## T2- rotation period of the second planet,
## a1 - large semi-axis of the orbit of the first planet,
## a2 - large semi-axis of the orbit of the second planet.

# TODO: remove law since a simpler one exists in `gravity` folder?

first_period = Symbol("first_period", units.time)
second_period = Symbol("second_period", units.time)
first_semi_axis = Symbol("first_semi_axis", units.length)
second_semi_axis = Symbol("second_semi_axis", units.length)

law = Eq((first_period)**2 / (second_period)**2, (first_semi_axis)**3 / (second_semi_axis)**3)


@validate_input(second_period_=second_period,
    first_semi_axis_=first_semi_axis,
    second_semi_axis_=second_semi_axis)
@validate_output(first_period)
def calculate_first_period(second_period_: Quantity, first_semi_axis_: Quantity,
    second_semi_axis_: Quantity) -> Quantity:
    result_expr = solve(law, first_period, dict=True)[1][first_period]
    result_expr = result_expr.subs({
        second_period: second_period_,
        first_semi_axis: first_semi_axis_,
        second_semi_axis: second_semi_axis_
    })
    return Quantity(result_expr)
