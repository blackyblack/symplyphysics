from sympy import (Eq, solve)
from symplyphysics import (
    symbols,
    units,
    Quantity,
    validate_input,
    validate_output,
    quantities,
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

lifetime = symbols.time
"""
Lifetime (:symbols:`time`) of the star on the main sequence.
"""

star_mass = symbols.mass
"""
:symbols:`mass` of the star.
"""

star_luminosity = symbols.luminocity
"""
:symbols:`luminocity` of the star.
"""

sun_lifetime = Quantity(1e10 * units.common_year, display_symbol="t_Sun", display_latex="t_\\odot")
"""

"""

sun_luminosity = Quantity(3.827e26 * units.watt, display_symbol="L_Sun", display_latex="L_\\odot")

law = Eq(lifetime,
    sun_lifetime * (star_mass / quantities.solar_mass) * (sun_luminosity / star_luminosity))


@validate_input(mass_of_star_=star_mass, luminosity_of_star_=star_luminosity)
@validate_output(lifetime)
def calculate_lifetime(mass_of_star_: Quantity, luminosity_of_star_: float) -> Quantity:
    result_expr = solve(law, lifetime, dict=True)[0][lifetime]
    result_expr = result_expr.subs({
        star_mass: mass_of_star_,
        star_luminosity: luminosity_of_star_,
    })
    return Quantity(result_expr)
