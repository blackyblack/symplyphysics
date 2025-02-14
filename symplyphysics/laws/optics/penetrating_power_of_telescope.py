from sympy import Eq, solve, log
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
    symbols,
)

# Description
## The penetrating power is the magnitude of the faintest stars visible with a telescope when observed at the zenith.

## Law is: m = 2.5 + 5 * lg(D), where
## m - penetrating power of telescope,
## D - diameter of the lens in millimeters.

# Possible link: https://ru.wikipedia.org/wiki/%D0%9F%D1%80%D0%BE%D0%BD%D0%B8%D1%86%D0%B0%D1%8E%D1%89%D0%B0%D1%8F_%D1%81%D0%B8%D0%BB%D0%B0
# TODO: find English link
# TODO: update documentation

penetrating_power_telescope = SymbolNew("penetrating_power_telescope", dimensionless)
# TODO: add to `symbols`

lens_diameter = symbols.diameter

one_millimeter = Quantity(1 * units.millimeter, display_symbol="D_0")

law = Eq(penetrating_power_telescope, 2.5 + 5 * log(lens_diameter / one_millimeter, 10))


@validate_input(lens_diameter_=lens_diameter)
@validate_output(penetrating_power_telescope)
def calculate_penetrating_power(lens_diameter_: Quantity) -> float:
    result_expr = solve(law, penetrating_power_telescope, dict=True)[0][penetrating_power_telescope]
    result_expr = result_expr.subs({
        lens_diameter: lens_diameter_,
    })
    return convert_to_float(result_expr)
