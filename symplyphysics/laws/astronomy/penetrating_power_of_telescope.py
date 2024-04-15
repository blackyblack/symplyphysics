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
## The penetrating force is the magnitude of the faintest stars visible with a telescope when observed at the zenith.

## Law is: m = 2.5 + 5 * lg(D), where
## m - penetrating power of telescope,
## D - diameter of the lens in millimeters.

penetrating_power_telescope = Symbol("penetrating_power_telescope", dimensionless)

lens_diameter = Symbol("lens_diameter", units.length)

one_millimeter = Quantity(1 * units.millimeter)

law = Eq(penetrating_power_telescope, 2.5 + 5 * log(lens_diameter / one_millimeter, 10))


def print_law() -> str:
    return print_expression(law)


@validate_input(lens_diameter_=lens_diameter)
@validate_output(penetrating_power_telescope)
def calculate_penetrating_power(lens_diameter_: Quantity) -> float:
    result_expr = solve(law, penetrating_power_telescope, dict=True)[0][penetrating_power_telescope]
    result_expr = result_expr.subs({
        lens_diameter: lens_diameter_,
    })
    return convert_to_float(result_expr)
