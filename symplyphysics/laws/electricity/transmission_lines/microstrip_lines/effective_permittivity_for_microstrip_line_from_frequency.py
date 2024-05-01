from sympy import Eq, solve, sqrt, log
from sympy.physics.units import speed_of_light
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float
)

## Description
## The microstrip line is a dielectric substrate on which a metal strip is applied.
## When a wave propagates along a microstrip line, part of the field goes out, since the microstrip line does
## not have metal borders on all sides, unlike, for example, rectangular waveguides. Then imagine an environment
## in which the field will have the same magnitude as the field of a microstrip line. The permittivity of such a
## medium will be called the effective permittivity of the line.

## Law is: eff = ((sqrt(er) - sqrt(ef)) / (1 + 4 * (Kf * f)^(-1.5)) + sqrt(ef))^2, where
## Kf = (4 * h * sqrt(er - 1) / c / 2) * (1 + 2 * ln(1 + W / h))^2,
## eff - effective permittivity of the microstrip line taking into account the dependence on frequency,
## ef - effective permittivity of the microstrip line without taking into account the dependence on frequency,
## er - relative permittivity of the dielectric substrate of the microstrip line,
## W - width of the microstrip line,
## h - thickness of substrate,
## f - frequency of signal,
## c - speed of light.

effective_permittivity = Symbol("effective_permittivity", dimensionless)

relative_permittivity = Symbol("relative_permittivity", dimensionless)
frequency = Symbol("frequency", units.frequency)
thickness_of_substrate = Symbol("thickness_of_substrate", units.length)
width = Symbol("width", units.length)
effective_permittivity_without_frequency = Symbol("effective_permittivity_without_frequency", dimensionless)

expression_1 = 4 * thickness_of_substrate * sqrt(relative_permittivity - 1) / speed_of_light / 2
expression_2 = (1 + 2 * log(1 + width / thickness_of_substrate))**2
expression_3 = expression_1 * expression_2
expression_4 = sqrt(relative_permittivity) - sqrt(effective_permittivity_without_frequency)

law = Eq(effective_permittivity, ((expression_4 / (1 + 4 * (expression_3 * frequency)**(-1.5))) + sqrt(effective_permittivity_without_frequency))**2)


@validate_input(relative_permittivity_=relative_permittivity,
    frequency_=frequency,
    thickness_of_substrate_=thickness_of_substrate,
    width_=width,
    effective_permittivity_without_frequency_=effective_permittivity_without_frequency)
@validate_output(effective_permittivity)
def calculate_effective_permittivity(relative_permittivity_: float,
    frequency_: Quantity, thickness_of_substrate_: Quantity,
    width_: Quantity, effective_permittivity_without_frequency_: float) -> float:
    result_expr = solve(law, effective_permittivity, dict=True)[0][effective_permittivity]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        frequency: frequency_,
        thickness_of_substrate: thickness_of_substrate_,
        width: width_,
        effective_permittivity_without_frequency: effective_permittivity_without_frequency_
    })
    return convert_to_float(result_expr)
