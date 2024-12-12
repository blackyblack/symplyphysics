from sympy import (Eq, solve)
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## The main sequence is the stage of stellar evolution.
## This stage begins in stars after the protostar stage. At the beginning of the main sequence stage, the age of the star
## is considered to be zero.
## It is possible to calculate the time spent on the main sequence for a star by knowing such a time for the Sun, as well
## as the mass and luminosity of the Sun and the star.

## Law is: T = ts * (M / Ms) * (Ls / L), where
## T - lifetime on the main sequence of the star,
## ts - lifetime on the main sequence of the Sun,
## M - mass of the star,
## Ms - mass of the Sun,
## Ls - luminosity of the Sun,
## L - luminosity of the star.

# Link: Wikipedia, second formula <https://en.wikipedia.org/wiki/Main_sequence#Lifetime>

lifetime = Symbol("lifetime", units.time)

mass_of_star = symbols.mass
luminosity_of_star = Symbol("luminosity_of_star", units.power)

mass_of_sun = Quantity(1.989e30 * units.kilogram)
lifetime_of_sun = Quantity(1e10 * units.common_year)
luminosity_of_sun = Quantity(3.827e26 * units.watt)

law = Eq(lifetime,
    lifetime_of_sun * (mass_of_star / mass_of_sun) * (luminosity_of_sun / luminosity_of_star))


@validate_input(mass_of_star_=mass_of_star, luminosity_of_star_=luminosity_of_star)
@validate_output(lifetime)
def calculate_lifetime(mass_of_star_: Quantity, luminosity_of_star_: float) -> Quantity:
    result_expr = solve(law, lifetime, dict=True)[0][lifetime]
    result_expr = result_expr.subs({
        mass_of_star: mass_of_star_,
        luminosity_of_star: luminosity_of_star_,
    })
    return Quantity(result_expr)
