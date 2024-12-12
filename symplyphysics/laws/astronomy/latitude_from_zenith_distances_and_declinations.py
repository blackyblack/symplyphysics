from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    angle_type,
)

# Description
## Knowing the zenith distances and declinations of the northern and southern stars, respectively, it is possible
## to determine the latitude of the observation site.
## Northern star is any star north of the zenith with known declination.
## Southern star is any star south of the zenith with known declination.

## Law is: phi = (zs - zn + ds + dn) / 2, where
## phi - latitude,
## zs - zenith distance for the north star,
## zn - zenith distance for the south star,
## ds - declination for the north star,
## dn - declination for the south star.

# Conditions:
# - both stars are at upper transit (culmination).

# TODO: find link

latitude = Symbol("zenith_distance_north", angle_type)

zenith_distance_north = Symbol("zenith_distance_north", angle_type)
zenith_distance_south = Symbol("zenith_distance_south", angle_type)
declination_north = Symbol("declination_north", angle_type)
declination_south = Symbol("declination_south", angle_type)

law = Eq(latitude,
    (zenith_distance_south - zenith_distance_north + declination_south + declination_north) / 2)


@validate_input(
    zenith_distance_north_=zenith_distance_north,
    zenith_distance_south_=zenith_distance_south,
    declination_north_=declination_north,
    declination_south_=declination_south,
)
@validate_output(latitude)
def calculate_latitude(zenith_distance_north_: Quantity, zenith_distance_south_: Quantity,
    declination_north_: Quantity, declination_south_: Quantity) -> Quantity:
    result_expr = solve(law, latitude, dict=True)[0][latitude]
    result_expr = result_expr.subs({
        zenith_distance_north: zenith_distance_north_,
        zenith_distance_south: zenith_distance_south_,
        declination_north: declination_north_,
        declination_south: declination_south_
    })
    return Quantity(result_expr)
