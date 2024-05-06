from sympy import Eq, solve, sqrt, pi, log
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
## The coplanar transmission line is a dielectric substrate on the surface of which 3 electrodes are located.
## The wave resistance of a transmission line is a value determined by the ratio of the voltage of the incident
## wave to the current of this wave in the transmission line.
## When a wave propagates along a coplanar line, part of the field goes out, since the coplanar line does
## not have metal borders on all sides, unlike, for example, rectangular waveguides. Then imagine an environment
## in which the field will have the same magnitude as the field of a coplanar line. The permittivity of such a
## medium will be called the effective permittivity of the line.

## Law is: Z = (30 * pi / sqrt(ef)) * (ln(2 * (1 + sqrt(sqrt(1 - k^2))) / (1 - sqrt(sqrt(1 - k^2)))) / pi), where
## k = a / b,
## Z - wave resistance of the coplanar line,
## ef - effective permittivity of the coplanar line,
## a - width of the central electrode of the coplanar line,
## b - distance between the extreme electrodes.

# Conditions:
# h - thickness of substrate.
# - h < b / 4;
# - 0 < k^2 <= 0.5.


wave_resistance = Symbol("wave_resistance", units.impedance)

effective_permittivity = Symbol("effective_permittivity", dimensionless)
distance_between_electrodes = Symbol("distance_between_electrodes", units.length)
central_electrode_width = Symbol("central_electrode_width", units.length)

constant_resistance = Quantity(30 * pi * units.ohm)
expression_1 = constant_resistance / sqrt(effective_permittivity)
expression_2 = sqrt(1 - (central_electrode_width / distance_between_electrodes)**2)
expression_3 =  log(2 * (1 + sqrt(expression_2)) / (1 - sqrt(expression_2))) / pi


law = Eq(wave_resistance, expression_1 * expression_3)


@validate_input(effective_permittivity_=effective_permittivity,
    distance_between_electrodes_=distance_between_electrodes,
    central_electrode_width_=central_electrode_width)
@validate_output(wave_resistance)
def calculate_wave_resistance(effective_permittivity_: float,
    distance_between_electrodes_: Quantity,
    central_electrode_width_: Quantity) -> float:
    if (central_electrode_width_.scale_factor / distance_between_electrodes_.scale_factor)**2 > 0.5:
        raise ValueError("k^2 must be less than or equal to the 0.5")

    result_expr = solve(law, wave_resistance, dict=True)[0][wave_resistance]
    result_expr = result_expr.subs({
        effective_permittivity: effective_permittivity_,
        distance_between_electrodes: distance_between_electrodes_,
        central_electrode_width: central_electrode_width_
    })
    return Quantity(result_expr)
