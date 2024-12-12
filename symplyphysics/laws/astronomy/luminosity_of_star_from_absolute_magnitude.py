from sympy import Eq, solve, log
from symplyphysics import (
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
)

# Description
## The luminosity of a star is related to the absolute magnitude of the star.
## The luminosity of a star is called the power of light energy emission compared to the power of light emission by the Sun.

## Law is: lg(L) = 0.4 * (5 - M), where
## L - luminosity of the star,
## M - absolute_magnitude of the star.

# Links: potential formula within text <https://en.wikipedia.org/wiki/Luminosity#Relationship_to_magnitude>

luminosity = Symbol("luminosity", dimensionless)

absolute_magnitude = Symbol("absolute_magnitude", dimensionless)

law = Eq(log(luminosity, 10), 0.4 * (5 - absolute_magnitude))


@validate_input(absolute_magnitude_=absolute_magnitude)
@validate_output(luminosity)
def calculate_luminosity(absolute_magnitude_: float) -> float:
    result_expr = solve(law, luminosity, dict=True)[0][luminosity]
    result_expr = result_expr.subs({
        absolute_magnitude: absolute_magnitude_,
    })
    return convert_to_float(result_expr)
