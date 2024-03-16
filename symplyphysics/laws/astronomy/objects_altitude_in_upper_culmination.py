from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    angle_type,
)

# Description
## In observational astronomy, culmination is the passage of a celestial object across the observer's local meridian.
## An object's altitude in degrees at its upper culmination is equal to 90 minus the observer's latitude plus
## the object's declination.
## https://учисьучись.рф/materials/shkolnaya-programma/astronomy/vysotasvetilvkulminacii/

## Law is: h = 90 - L + d, where
## h - object's altitude,
## L - latitude,
## d - object's declination.

altitude = Symbol("altitude", angle_type)

latitude = Symbol("latitude", angle_type)
declination = Symbol("declination", angle_type)

angle_constant = Quantity(90 * units.deg)

law = Eq(altitude, angle_constant - latitude + declination)


def print_law() -> str:
    return print_expression(law)


@validate_input(latitude_=latitude, declination_=declination)
@validate_output(altitude)
def calculate_altitude(latitude_: Quantity, declination_: Quantity) -> Quantity:
    result_expr = solve(law, altitude, dict=True)[0][altitude]
    result_expr = result_expr.subs({
        latitude: latitude_,
        declination: declination_,
    })
    return Quantity(result_expr)
