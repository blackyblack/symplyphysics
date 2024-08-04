from sympy import Eq, solve, log, pi
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
## The effective width of a microstrip line is the width of such a flat capacitor, the electric intensity between the plates
## of which is equal to the electric intensity in the dielectric of the substrate under the line strip.

## Law is: Wef / h = W / h + (1.25 * t / (pi * h)) * (1 + ln(2 * h / t)), where
## Wef - effective width of the microstrip line,
## W - width of the microstrip line,
## h - thickness of substrate,
## t - strip thickness of the microstrip line.

# Conditions:
# - the thickness of the substrate of the microstrip line should be less than the width * 2 * pi.

effective_width = Symbol("effective_width", units.length)

strip_thickness = Symbol("strip_thickness", units.length)
thickness_of_substrate = Symbol("thickness_of_substrate", units.length)
width = Symbol("width", units.length)

expression_1 = width / thickness_of_substrate
expression_2 = 1.25 * strip_thickness / (pi * thickness_of_substrate)
expression_3 = 1 + log(2 * thickness_of_substrate / strip_thickness)

law = Eq(effective_width / thickness_of_substrate, expression_1 + expression_2 * expression_3)


def print_law() -> str:
    return print_expression(law)


@validate_input(strip_thickness_=strip_thickness,
    thickness_of_substrate_=thickness_of_substrate,
    width_=width)
@validate_output(effective_width)
def calculate_effective_width(strip_thickness_: Quantity, thickness_of_substrate_: Quantity,
    width_: Quantity) -> Quantity:
    if thickness_of_substrate_.scale_factor >= width_.scale_factor * 2 * pi:
        raise ValueError("The thickness of substrate must be less than the width * 2 * pi")
    result_expr = solve(law, effective_width, dict=True)[0][effective_width]
    result_expr = result_expr.subs({
        strip_thickness: strip_thickness_,
        thickness_of_substrate: thickness_of_substrate_,
        width: width_
    })
    return Quantity(result_expr)
