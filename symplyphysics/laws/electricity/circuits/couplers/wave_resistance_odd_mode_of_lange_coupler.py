from sympy import Eq, solve, sqrt
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
)

## Description
## The Lange coupler is based on microstrip transmission lines. When this coupler is in operation, both even
## and odd modes are distributed. Knowing the coupling coefficient between the coupler segments, the characteristic
## resistance of the transmission line to which the coupler is connected, as well as the number of coupler segments,
## it is possible to calculate the wave resistance for an odd mode.
## https://habrastorage.org/r/w1560/getpro/habr/upload_files/054/d02/c8d/054d02c8d91c06425ae079d34b18ce15.jpeg

## Law is: Z = (R0 * sqrt((1 - C0) / (1 + C0))) * (N - 1) * (1 + sqrt(C0**2 + (1 - C0**2) * (N - 1)**2)) / (C0 + sqrt(C0**2 + (1 - C0**2) * (N - 1)**2) + (N - 1) * (1 - C0)), where
## Z - wave resistance of the odd mode,
## R0 - characteristic resistance of the transmission line,
## C0 - coupling factor between the coupler segments,
## N - the number of segments of the Lange coupler.

wave_resistance_odd_modes = Symbol("wave_resistance_odd_modes", units.impedance)

coupling_factor = Symbol("coupling_factor", dimensionless)
characteristic_resistance = Symbol("characteristic_resistance", units.impedance)
number_segments = Symbol("number_segments", dimensionless)

expression_1 = characteristic_resistance * sqrt((1 - coupling_factor) / (1 + coupling_factor))
expression_2 = sqrt(coupling_factor**2 + (1 - coupling_factor**2) * (number_segments - 1)**2)
expression_3 = (number_segments - 1) * (1 + expression_2) / (coupling_factor + expression_2 +
    (number_segments - 1) * (1 - coupling_factor))
law = Eq(wave_resistance_odd_modes, expression_1 * expression_3)


def print_law() -> str:
    return print_expression(law)


@validate_input(coupling_factor_=coupling_factor,
    characteristic_resistance_=characteristic_resistance,
    number_segments_=number_segments)
@validate_output(wave_resistance_odd_modes)
def calculate_wave_resistance_odd_modes(coupling_factor_: float,
    characteristic_resistance_: Quantity, number_segments_: int) -> Quantity:
    result_expr = solve(law, wave_resistance_odd_modes, dict=True)[0][wave_resistance_odd_modes]
    result_expr = result_expr.subs({
        coupling_factor: coupling_factor_,
        characteristic_resistance: characteristic_resistance_,
        number_segments: number_segments_
    })
    return Quantity(result_expr)
