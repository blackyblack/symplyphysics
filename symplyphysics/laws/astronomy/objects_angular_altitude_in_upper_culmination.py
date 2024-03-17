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
## An object's angular altitude in degrees at its upper culmination is equal to 90 minus the observer's latitude plus
## the object's declination.
## Object's declination is equal to the angular distance on the celestial sphere from the plane of the celestial equator to the sun.
## https://учисьучись.рф/materials/shkolnaya-programma/astronomy/vysotasvetilvkulminacii/

## Law is: h = 90 - L + d, where
## h - object's angular altitude,
## L - latitude,
## d - object's declination.

# Conditions:
# - it is valid when object is at its upper culmination.

angular_altitude = Symbol("angular_altitude", angle_type)

latitude = Symbol("latitude", angle_type)
declination = Symbol("declination", angle_type)

ninety_degrees = Quantity(90 * units.deg)

law = Eq(angular_altitude, ninety_degrees - latitude + declination)


def print_law() -> str:
    return print_expression(law)


@validate_input(latitude_=latitude, declination_=declination)
@validate_output(angular_altitude)
def calculate_angular_altitude(latitude_: Quantity, declination_: Quantity) -> Quantity:
    result_expr = solve(law, angular_altitude, dict=True)[0][angular_altitude]
    result_expr = result_expr.subs({
        latitude: latitude_,
        declination: declination_,
    })
    return Quantity(result_expr)
