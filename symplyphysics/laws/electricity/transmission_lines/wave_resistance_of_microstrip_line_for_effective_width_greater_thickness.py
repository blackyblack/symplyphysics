from sympy import Eq, solve, sqrt, log, pi
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless
)

## Description
## The microstrip line is a dielectric substrate on which a metal strip is applied.
## The wave resistance of a transmission line is a value determined by the ratio of the voltage of the incident
## wave to the current of this wave in the transmission line.
## When a wave propagates along a microstrip line, part of the field goes out, since the microstrip line does
## not have metal borders on all sides, unlike, for example, rectangular waveguides. Then imagine an environment
## in which the field will have the same magnitude as the field of a microstrip line. The permittivity of such a
## medium will be called the effective permittivity of the line.
## The effective width of a microstrip line is the width of such a flat capacitor, the electric intensity between the plates
## of which is equal to the electric intensity in the dielectric of the substrate under the line strip.

## Law is: Z = 120 * pi * (Wef / h + 1.393 + 0.667 * ln(Wef / h + 1.444))**(-1) / sqrt(ef), where
## Z - wave resistance of the microstrip line,
## Wef - effective width of the microstrip line,
## h - thickness of substrate,
## ef - effective permittivity of the microstrip line.

# Conditions:
# - the effective width of the microstrip line should be greater than the thickness of the substrate.

resistance = Symbol("resistance", units.impedance)

effective_permittivity = Symbol("effective_permittivity", dimensionless)
thickness_of_substrate = Symbol("thickness_of_substrate", units.length)
effective_width = Symbol("effective_width", units.length)

expression_1 = effective_width / thickness_of_substrate + 1.393
expression_2 = 0.667 * log(effective_width / thickness_of_substrate + 1.444)
expression_3 = (expression_1 + expression_2)**(-1)
constant_resistance = Quantity(120 * pi * units.ohm)
law = Eq(resistance, constant_resistance * expression_3 / sqrt(effective_permittivity))


def print_law() -> str:
    return print_expression(law)


@validate_input(effective_permittivity_=effective_permittivity,
    thickness_of_substrate_=thickness_of_substrate,
    effective_width_=effective_width)
@validate_output(resistance)
def calculate_resistance(effective_permittivity_: float, thickness_of_substrate_: Quantity,
    effective_width_: Quantity) -> Quantity:
    if thickness_of_substrate_.scale_factor > effective_width_.scale_factor:
        raise ValueError("The thickness of substrate must be less than the effective width")
    result_expr = solve(law, resistance, dict=True)[0][resistance]
    result_expr = result_expr.subs({
        effective_permittivity: effective_permittivity_,
        thickness_of_substrate: thickness_of_substrate_,
        effective_width: effective_width_
    })
    return Quantity(result_expr)
