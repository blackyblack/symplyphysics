from sympy import Eq, solve, pi, log
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

## Description
## The microstrip line is a dielectric substrate on which a metal strip is applied.
## The effective width of a microstrip line is the width of such a flat capacitor, the electric intensity between the plates
## of which is equal to the electric intensity in the dielectric of the substrate under the line strip.
## The attenuation coefficient shows how many times the transmitted signal weakens per unit length of the microstrip line.

## Law is: am = (1.38 * Rs / (h * Z0)) * ((32 - (Wef / h)^2) / (32 + (Wef / h)^2)) * F, where
## F = 1 + (h / Wef) * (1 - 1.25 * t / (pi * h) + 1.25 * ln(2 * h / t) / pi) (if h < 2 * pi * Wef),
## am - attenuation coefficient of the metal of the microstrip line,
## Wef - effective width of the microstrip line,
## W - width of the microstrip line,
## Rs - surface resistance of the metal strip,
## h - thickness of substrate,
## Z0 - wave resistance of the microstrip line,
## t - strip thickness of the microstrip line.

# Conditions:
# - the thickness of the substrate of the microstrip line should be greater than or equal to the effective width;
# - the thickness of the substrate of the microstrip line should be less than the effective width * 2 * pi.

attenuation_coefficient = Symbol("attenuation_coefficient", 1 / units.length)

surface_resistance = Symbol("surface_resistance", units.impedance)
wave_resistance = Symbol("wave_resistance", units.impedance)
thickness_of_substrate = Symbol("thickness_of_substrate", units.length)
effective_width = Symbol("effective_width", units.length)
strip_thickness = Symbol("strip_thickness", units.length)
width = Symbol("width", units.length)

expression_1 = (effective_width / thickness_of_substrate)**2
expression_2 = 1.38 * surface_resistance / (thickness_of_substrate * wave_resistance)
expression_3 = 1 + (thickness_of_substrate / effective_width) * (1 - 1.25 * strip_thickness /
    (pi * thickness_of_substrate) + 1.25 * log(2 * thickness_of_substrate / strip_thickness) / pi)

law = Eq(attenuation_coefficient,
    expression_2 * ((32 - expression_1) / (32 + expression_1)) * expression_3)


@validate_input(surface_resistance_=surface_resistance,
    wave_resistance_=wave_resistance,
    thickness_of_substrate_=thickness_of_substrate,
    effective_width_=effective_width,
    strip_thickness_=strip_thickness,
    width_=width)
@validate_output(attenuation_coefficient)
def calculate_attenuation_coefficient(surface_resistance_: Quantity, wave_resistance_: Quantity,
    thickness_of_substrate_: Quantity, effective_width_: Quantity, strip_thickness_: Quantity,
    width_: Quantity) -> Quantity:
    if thickness_of_substrate_.scale_factor < effective_width_.scale_factor:
        raise ValueError(
            "The thickness of substrate must be greater than or equal to the effective width")
    if thickness_of_substrate_.scale_factor >= effective_width_.scale_factor * 2 * pi:
        raise ValueError(
            "The thickness of substrate must be less than the effective width * 2 * pi")

    result_expr = solve(law, attenuation_coefficient, dict=True)[0][attenuation_coefficient]
    result_expr = result_expr.subs({
        surface_resistance: surface_resistance_,
        wave_resistance: wave_resistance_,
        thickness_of_substrate: thickness_of_substrate_,
        effective_width: effective_width_,
        strip_thickness: strip_thickness_,
        width: width_,
    })
    return Quantity(result_expr)
