from sympy import Eq, solve, log
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

## Description
## The microstrip line is a dielectric substrate on which a metal strip is applied.

## Law is: R = (1.4 + 0.217 * ln(W / (5 * t))) * Rs * l / (2 * (W + t)), where
## R - resistance of the microstrip line,
## W - strip width of the microstrip line,
## l - strip length,
## t - strip thickness,
## Rs - surface resistance of the metal strip.

resistance = Symbol("resistance", units.impedance)

strip_thickness = Symbol("strip_thickness", units.length)
strip_length = Symbol("strip_length", units.length)
strip_width = Symbol("strip_width", units.length)
surface_resistance = Symbol("surface_resistance", units.impedance)

expression_1 = surface_resistance * strip_length / (2 * (strip_width + strip_thickness))

law = Eq(resistance, (1.4 + 0.217 * log(strip_width / (5 * strip_thickness))) * expression_1)


def print_law() -> str:
    return print_expression(law)


@validate_input(strip_thickness_=strip_thickness,
    strip_length_=strip_length,
    strip_width_=strip_width,
    surface_resistance_=surface_resistance)
@validate_output(resistance)
def calculate_resistance(strip_thickness_: Quantity, strip_length_: Quantity,
    strip_width_: Quantity, surface_resistance_: Quantity) -> Quantity:
    result_expr = solve(law, resistance, dict=True)[0][resistance]
    result_expr = result_expr.subs({
        strip_thickness: strip_thickness_,
        strip_length: strip_length_,
        strip_width: strip_width_,
        surface_resistance: surface_resistance_,
    })
    return Quantity(result_expr)
