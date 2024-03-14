from sympy import (Eq, solve, S)
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
## The main sequence is the stage of stellar evolution formed by the stars at this stage and the corresponding luminosity class.
## The main sequence of stars falls after the protostar stage. At this point, the age of the star is considered zero and it is
## located on the so-called initial main sequence.

## Law is: T = 10 * M^(1 - n), where
## T - lifetime on the main sequence in billions of years,
## M - mass of the star (in masses of the Sun),
## n - indicator depending on the mass of the star.

lifetime = Symbol("lifetime", dimensionless)

mass_of_star = Symbol("mass_of_star", dimensionless)
indicator = Symbol("indicator", dimensionless)

law = Eq(lifetime, 10 * mass_of_star**(1 - indicator))


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_of_star_=mass_of_star, indicator_=indicator)
@validate_output(lifetime)
def calculate_lifetime(mass_of_star_: float, indicator_: float) -> float:
    result_expr = solve(law, lifetime, dict=True)[0][lifetime]
    result_expr = result_expr.subs({
        mass_of_star: mass_of_star_,
        indicator: indicator_,
    })
    return float(convert_to(Quantity(result_expr), S.One).evalf())
