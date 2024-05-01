from sympy import Eq, solve, sqrt
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

## Law is: Zf = Z0 * sqrt(ef / eff) * (eff - 1) / (ef - 1), where
## Zf - wave resistance of the microstrip line taking into account the dependence on frequency,
## Z0 - wave resistance of the microstrip line without taking into account the dependence on frequency,
## eff - effective permittivity of the microstrip line taking into account the dependence on frequency,
## ef - effective permittivity of the microstrip line without taking into account the dependence on frequency.

wave_resistance = Symbol("wave_resistance", units.impedance)

wave_resistance_without_frequency = Symbol("wave_resistance_without_frequency", units.impedance)
effective_permittivity = Symbol("effective_permittivity", dimensionless)
effective_permittivity_without_frequency = Symbol("effective_permittivity_without_frequency", dimensionless)

expression_1 = sqrt(effective_permittivity_without_frequency / effective_permittivity)
expression_2 = (effective_permittivity - 1) / (effective_permittivity_without_frequency - 1)

law = Eq(wave_resistance, wave_resistance_without_frequency * expression_1 * expression_2)


def print_law() -> str:
    return print_expression(law)


@validate_input(wave_resistance_without_frequency_=wave_resistance_without_frequency,
    effective_permittivity_=effective_permittivity,
    effective_permittivity_without_frequency_=effective_permittivity_without_frequency)
@validate_output(wave_resistance)
def calculate_wave_resistance(wave_resistance_without_frequency_: Quantity, effective_permittivity_: float,
    effective_permittivity_without_frequency_: float) -> Quantity:
    result_expr = solve(law, wave_resistance, dict=True)[0][wave_resistance]
    result_expr = result_expr.subs({
        wave_resistance_without_frequency: wave_resistance_without_frequency_,
        effective_permittivity: effective_permittivity_,
        effective_permittivity_without_frequency: effective_permittivity_without_frequency_
    })
    return Quantity(result_expr)
