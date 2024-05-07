from sympy import Eq, solve, pi, log
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless
)

## Description
## The microstrip line is a dielectric substrate on which a metal strip is applied.
## The effective width of a microstrip line is the width of such a flat capacitor, the electric intensity between the plates
## of which is equal to the electric intensity in the dielectric of the substrate under the line strip.
## The attenuation coefficient shows how many times the transmitted signal weakens per unit length of the microstrip line.
## When a wave propagates along a microstrip line, part of the field goes out, since the microstrip line does
## not have metal borders on all sides, unlike, for example, rectangular waveguides. Then imagine an environment
## in which the field will have the same magnitude as the field of a microstrip line. The permittivity of such a
## medium will be called the effective permittivity of the line.

## Law is: am = (6.1e-5 * Rs * Z0 * ef / h) * (Wef / h + (0.667 * Wef / h) / (Wef /  h + 1.444)) * F, where
## F = 1 + (h / Wef) * (1 - 1.25 * t / (pi * h) + 1.25 * ln(2 * h / t) / pi),
## am - attenuation coefficient of the metal of the microstrip line,
## Wef - effective width of the microstrip line,
## Rs - surface resistance of the metal strip,
## h - thickness of substrate,
## Z0 - wave resistance of the microstrip line,
## t - strip thickness of the microstrip line,
## ef - effective permittivity of the microstrip line.

# Conditions:
# - the thickness of the substrate of the microstrip line should be less than the effective width.

attenuation_coefficient = Symbol("attenuation_coefficient", 1 / units.length)

surface_resistance = Symbol("surface_resistance", units.impedance)
wave_resistance = Symbol("wave_resistance", units.impedance)
thickness_of_substrate = Symbol("thickness_of_substrate", units.length)
effective_width = Symbol("effective_width", units.length)
strip_thickness = Symbol("strip_thickness", units.length)
effective_permittivity = Symbol("effective_permittivity", dimensionless)

expression_1 = effective_width / thickness_of_substrate
constant = Quantity(6.1e-5 * units.ohm**(-2))
expression_2 = constant * surface_resistance * wave_resistance * effective_permittivity / thickness_of_substrate
expression_3 = 1 + (thickness_of_substrate / effective_width) * (1 - 1.25 * strip_thickness / (pi * thickness_of_substrate) + 1.25 * log(2 * thickness_of_substrate / strip_thickness) / pi)


law = Eq(attenuation_coefficient, expression_2 * (expression_1 + 0.667 * expression_1 / (expression_1 + 1.444)) * expression_3)


@validate_input(surface_resistance_=surface_resistance,
    wave_resistance_=wave_resistance,
    thickness_of_substrate_=thickness_of_substrate,
    effective_width_=effective_width,
    strip_thickness_=strip_thickness,
    effective_permittivity_=effective_permittivity)
@validate_output(attenuation_coefficient)
def calculate_attenuation_coefficient(surface_resistance_: Quantity,
    wave_resistance_: Quantity, thickness_of_substrate_: Quantity,
    effective_width_: Quantity, strip_thickness_: Quantity,
    effective_permittivity_: float) -> Quantity:
    if thickness_of_substrate_.scale_factor >= effective_width_.scale_factor:
        raise ValueError("The thickness of substrate must be less than the effective width")

    result_expr = solve(law, attenuation_coefficient, dict=True)[0][attenuation_coefficient]
    result_expr = result_expr.subs({
        surface_resistance: surface_resistance_,
        wave_resistance: wave_resistance_,
        thickness_of_substrate: thickness_of_substrate_,
        effective_width: effective_width_,
        strip_thickness: strip_thickness_,
        effective_permittivity: effective_permittivity_,
    })
    return Quantity(result_expr)
