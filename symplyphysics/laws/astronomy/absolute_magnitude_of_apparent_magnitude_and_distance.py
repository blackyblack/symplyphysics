from sympy import Eq, solve, log
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
)

# Description
## The absolute magnitude for stars is defined as the apparent magnitude of an object if it were located at a distance
## of 10 parsecs (2.063e+6 astronomical units) from the observer and would not experience either interstellar or atmospheric absorption.
## The apparent magnitude is a measure of the brightness of a celestial body (more precisely, the illumination created
## by this body) from the observer's point of view. The brighter the object, the smaller its magnitude.

## Law is: M = m - 5 * lg(d / d0), where
## M - absolute magnitude,
## m - apparent magnitude,
## d - distance to the object,
## d0 - constant equal to 2.063e+6 astronomical units.

absolute_magnitude = Symbol("absolute_magnitude", dimensionless)

apparent_magnitude = Symbol("apparent_magnitude", dimensionless)
distance = Symbol("distance", units.length)

distance_constant = Quantity(2.063e+6 * units.astronomical_unit)

law = Eq(absolute_magnitude, apparent_magnitude - 5 * log(distance / distance_constant, 10))


def print_law() -> str:
    return print_expression(law)


@validate_input(apparent_magnitude_=apparent_magnitude, distance_=distance)
@validate_output(absolute_magnitude)
def calculate_absolute_magnitude(apparent_magnitude_: float, distance_: Quantity) -> float:
    result_expr = solve(law, absolute_magnitude, dict=True)[0][absolute_magnitude]
    result_expr = result_expr.subs({
        apparent_magnitude: apparent_magnitude_,
        distance: distance_,
    })
    return convert_to_float(Quantity(result_expr))
