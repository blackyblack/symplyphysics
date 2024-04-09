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
## The apparent magnitude is a measure of the brightness of a celestial body (more precisely, the illumination created
## by this body) from the observer's point of view. The brighter the object, the smaller its magnitude.
## The relationship of the stellar magnitude scale with real physical quantities is logarithmic, since a change in brightness
## by the same number of times is perceived by the eye as a change by the same amount.
## The difference in the stellar magnitudes of two objects is equal to the decimal logarithm of the ratio of their illuminances,
## up to a multiplier.

## Law is: m2 - m1 = -2.5 * lg(L2 / L1), where
## m2 - apparent magnitude of second object,
## m1 - apparent magnitude of first object,
## l2 - illuminance of second object,
## l1 - illuminance of first object.

apparent_magnitude_first = Symbol("apparent_magnitude_first", dimensionless)
apparent_magnitude_second = Symbol("apparent_magnitude_second", dimensionless)
illuminance_first = Symbol("illuminance_first", units.energy / units.area)
illuminance_second = Symbol("illuminance_second", units.energy / units.area)

law = Eq(apparent_magnitude_second - apparent_magnitude_first,
    -2.5 * log(illuminance_second / illuminance_first, 10))


def print_law() -> str:
    return print_expression(law)


@validate_input(apparent_magnitude_first_=apparent_magnitude_first,
    illuminance_first_=illuminance_first,
    illuminance_second_=illuminance_second)
@validate_output(apparent_magnitude_second)
def calculate_apparent_magnitude_second(apparent_magnitude_first_: float,
    illuminance_first_: Quantity, illuminance_second_: Quantity) -> float:
    result_expr = solve(law, apparent_magnitude_second, dict=True)[0][apparent_magnitude_second]
    result_expr = result_expr.subs({
        apparent_magnitude_first: apparent_magnitude_first_,
        illuminance_first: illuminance_first_,
        illuminance_second: illuminance_second_
    })
    return convert_to_float(Quantity(result_expr))
