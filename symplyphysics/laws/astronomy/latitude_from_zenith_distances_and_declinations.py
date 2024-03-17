from sympy import (Eq, solve, log, S)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    convert_to,
    angle_type,
)

# Description
## Knowing the zenith distances and declinations of the northern and southern stars, respectively, it is possible
## to determine the latitude of the observation site.

## Law is: phi = (zs - zn + ds + dn) / 2, where
## phi - latitude,
## zs - zenith distance for the north star,
## zn - zenith distance for the south star,
## ds - declination for the north star,
## dn - declination for the south star.

latitude = Symbol("zenith_distance_north", angle_type)

zenith_distance_north = Symbol("zenith_distance_north", angle_type)
zenith_distance_south = Symbol("zenith_distance_south", angle_type)
declination_north = Symbol("declination_north", angle_type)
declination_south = Symbol("declination_south", angle_type)

law = Eq(latitude, (zenith_distance_south - zenith_distance_north + declination_south + declination_north) / 2)


def print_law() -> str:
    return print_expression(law)


@validate_input(zenith_distance_north_=zenith_distance_north,
    zenith_distance_south_=zenith_distance_south,
    declination_north_=declination_north,
    declination_south_=declination_south,)
@validate_output(latitude)
def calculate_latitude(zenith_distance_north_: Quantity, zenith_distance_south_: Quantity, declination_north_: Quantity,
    declination_south_: Quantity) -> Quantity:
    result_expr = solve(law, latitude, dict=True)[0][latitude]
    result_expr = result_expr.subs({
        zenith_distance_north: zenith_distance_north_,
        zenith_distance_south: zenith_distance_south_,
        declination_north: declination_north_,
        declination_south: declination_south_
    })
    return Quantity(result_expr)
