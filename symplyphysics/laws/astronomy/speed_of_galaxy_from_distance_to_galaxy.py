from sympy import (
    Eq,
    solve,
)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## Hubble's Law is a cosmological law describing the expansion of the universe.  It says that the speed of removal of a galaxy
## depends on the distance to it.

## Law is: v = H * r, where
## v - speed of galaxy (it is directed away from the observer looking at the galaxy),
## H - Hubble's constant,
## r - distance to galaxy.

# Links: Wikipedia <https://en.wikipedia.org/wiki/Hubble%27s_law#>

speed_of_galaxy = Symbol("speed_of_galaxy", units.velocity)

distance_to_galaxy = Symbol("distance_to_galaxy", units.length)

hubble_constant = Quantity(2.2e-18 / units.second)

law = Eq(speed_of_galaxy, hubble_constant * distance_to_galaxy)


@validate_input(distance_to_galaxy_=distance_to_galaxy)
@validate_output(speed_of_galaxy)
def calculate_speed(distance_to_galaxy_: Quantity) -> Quantity:
    result_expr = solve(law, speed_of_galaxy, dict=True)[0][speed_of_galaxy]
    result_expr = result_expr.subs({
        distance_to_galaxy: distance_to_galaxy_,
    })
    return Quantity(result_expr)
