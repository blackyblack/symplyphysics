from sympy import Eq, solve, sqrt
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, dimensionless,
    convert_to_float)

## Description
## The microstrip line is a dielectric substrate on which a metal strip is applied.
## When a wave propagates along a microstrip line, part of the field goes out, since the microstrip line does
## not have metal borders on all sides, unlike, for example, rectangular waveguides. Then imagine an environment
## in which the field will have the same magnitude as the field of a microstrip line. The permittivity of such a
## medium will be called the effective permittivity of the line.

## Law is: ef = (1 + er) / 2 + ((er - 1) / 2) * (1 + 12 * h / W)^(-1 / 2) - ((er - 1) / 4.6) * (t / h) * sqrt(h / W), where
## ef - effective permittivity of the microstrip line,
## er - relative permittivity of the dielectric substrate of the microstrip line,
## W - width of the microstrip line,
## h - thickness of substrate,
## t - strip thickness of the microstrip line.

# Conditions:
# - the thickness of the substrate of the microstrip line should be less than the width.

effective_permittivity = Symbol("effective_permittivity", dimensionless)

relative_permittivity = Symbol("relative_permittivity", dimensionless)
strip_thickness = Symbol("strip_thickness", units.length)
thickness_of_substrate = Symbol("thickness_of_substrate", units.length)
width = Symbol("width", units.length)

expression_1 = (1 + relative_permittivity) / 2
expression_2 = (relative_permittivity - 1) / 2
expression_3 = (1 + 12 * thickness_of_substrate / width)**(-1 / 2)
expression_4 = ((relative_permittivity - 1) / 4.6) * (strip_thickness /
    thickness_of_substrate) * sqrt(thickness_of_substrate / width)

law = Eq(effective_permittivity, expression_1 + expression_2 * expression_3 - expression_4)


@validate_input(relative_permittivity_=relative_permittivity,
    strip_thickness_=strip_thickness,
    thickness_of_substrate_=thickness_of_substrate,
    width_=width)
@validate_output(effective_permittivity)
def calculate_effective_permittivity(relative_permittivity_: float, strip_thickness_: Quantity,
    thickness_of_substrate_: Quantity, width_: Quantity) -> float:
    if thickness_of_substrate_.scale_factor >= width_.scale_factor:
        raise ValueError("The thickness of substrate must be less than the width")
    result_expr = solve(law, effective_permittivity, dict=True)[0][effective_permittivity]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        strip_thickness: strip_thickness_,
        thickness_of_substrate: thickness_of_substrate_,
        width: width_
    })
    return convert_to_float(result_expr)
