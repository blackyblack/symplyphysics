from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Kepler's third law relates the ratio of the rotation periods of the planets to the ratio of the large semi-axes
## of the orbits. Newton also clarified this law by adding the mass of the planet and the mass of the Sun.

## Law is: ((M + m1) * T1^2) / ((M + m2) * T2^2) = a1^3 / a2^3, where
## T1 - rotation period of the first planet,
## T2- rotation period of the second planet,
## a1 - large semi-axis of the orbit of the first planet,
## a2 - large semi-axis of the orbit of the second planet,
## M - mass of the Sun,
## m1 - mass of the first planet,
## m2 - mass of the second planet.

first_period = Symbol("first_period", units.time)
second_period = Symbol("second_period", units.time)
first_semi_axis = Symbol("first_semi_axis", units.length)
second_semi_axis = Symbol("second_semi_axis", units.length)
first_mass = Symbol("first_mass", units.mass)
second_mass = Symbol("second_mass", units.mass)

mass_of_sun = Quantity(1.989e30 * units.kilogram)

law = Eq(((mass_of_sun + first_mass) * first_period**2) / ((mass_of_sun + second_mass) * second_period**2), first_semi_axis**3 / second_semi_axis**3)


def print_law() -> str:
    return print_expression(law)


@validate_input(second_period_=second_period,
    first_semi_axis_=first_semi_axis,
    second_semi_axis_=second_semi_axis,
    first_mass_=first_mass, second_mass_=second_mass)
@validate_output(first_period)
def calculate_first_period(second_period_: Quantity, first_semi_axis_: Quantity,
    second_semi_axis_: Quantity, first_mass_: Quantity, second_mass_: Quantity) -> Quantity:
    result_expr = solve(law, first_period, dict=True)[1][first_period]
    result_expr = result_expr.subs({
        second_period: second_period_,
        first_semi_axis: first_semi_axis_,
        second_semi_axis: second_semi_axis_,
        first_mass: first_mass_,
        second_mass: second_mass_,
    })
    return Quantity(result_expr)
