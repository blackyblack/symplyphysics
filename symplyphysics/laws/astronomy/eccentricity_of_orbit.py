from sympy import Eq, solve, sqrt
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
)
from symplyphysics.core.convert import convert_to_dimensionless

# Description
## Many bodies in space move in elliptical orbits. For an elliptical orbit, the eccentricity can be calculated.
## The eccentricity can be used to judge the elongation of the elliptical orbit.

## Law is: e = sqrt(1 - (b / a)^2), where
## e - eccentricity  of orbit,
## b - small semi-axis of ellipse,
## a - large semi-axis of ellipse.

eccentricity = Symbol("eccentricity", dimensionless)

small_semi_axis = Symbol("small_semi_axis", units.length)
large_semi_axis = Symbol("large_semi_axis", units.length)

law = Eq(eccentricity, sqrt(1 - (small_semi_axis / large_semi_axis)**2))


def print_law() -> str:
    return print_expression(law)


@validate_input(small_semi_axis_=small_semi_axis, large_semi_axis_=large_semi_axis)
@validate_output(eccentricity)
def calculate_eccentricity(small_semi_axis_: Quantity, large_semi_axis_: Quantity) -> float:
    if small_semi_axis_.scale_factor > large_semi_axis_.scale_factor:
        raise ValueError("The small semi-axis must be less than or equal to the large semi-axis")
    result_expr = solve(law, eccentricity, dict=True)[0][eccentricity]
    result_expr = result_expr.subs({
        small_semi_axis: small_semi_axis_,
        large_semi_axis: large_semi_axis_,
    })
    return convert_to_dimensionless(Quantity(result_expr))
