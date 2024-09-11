from sympy import Eq, solve
from symplyphysics import (
    clone_symbol,
    symbols,
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
)

# Description
## The main sequence is the stage of stellar evolution.
## This stage begins in stars after the protostar stage. At the beginning of the main sequence stage, the age of the star
## is considered to be zero.
## The indicator is 4.75 for stars with a mass of 0.7 - 2 solar masses and is equal to 4.75 * M + 2.125 for stars with a mass
## of 0.1 – 0.7 solar masses.

## Law is: T = 10 * M^(1 - n), where
## T - lifetime on the main sequence in billions of years,
## M - mass of the star (in masses of the Sun),
## n - indicator depending on the mass of the star.

lifetime = Symbol("lifetime", units.time)
mass_of_star = clone_symbol(symbols.mass)
indicator = Symbol("indicator", dimensionless)

mass_of_sun = Quantity(1.989e30 * units.kilogram)
one_billion_years = Quantity(1e9 * units.common_year)

law = Eq(lifetime, 10 * one_billion_years * ((mass_of_star / mass_of_sun)**(1 - indicator)))


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_of_star_=mass_of_star, indicator_=indicator)
@validate_output(lifetime)
def calculate_lifetime(mass_of_star_: Quantity, indicator_: float) -> Quantity:
    result_expr = solve(law, lifetime, dict=True)[0][lifetime]
    result_expr = result_expr.subs({
        mass_of_star: mass_of_star_,
        indicator: indicator_,
    })
    return Quantity(result_expr)
