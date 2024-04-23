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
## and odd modes are distributed. Knowing the wave resistance for even and odd modes, it is possible to calculate
## the equivalent wave resistance of the coupler.
## https://habrastorage.org/r/w1560/getpro/habr/upload_files/054/d02/c8d/054d02c8d91c06425ae079d34b18ce15.jpeg

## Law is: Z = sqrt((Zo * Ze * (Zo + Ze)^2) / ((Zo + Ze * (N - 1)) * (Ze + Zo * (N - 1)))), where
## Z - equivalent wave resistance of the coupler,
## Zo - wave resistance of the odd mode,
## Ze - wave resistance of the even mode,
## N - the number of segments of the Lange coupler.

wave_resistance = Symbol("wave_resistance", units.impedance)

wave_resistance_odd_modes = Symbol("wave_resistance_odd_modes", units.impedance)
wave_resistance_even_modes = Symbol("wave_resistance_even_modes", units.impedance)
number_segments = Symbol("number_segments", dimensionless)

law = Eq(
    wave_resistance,
    sqrt((wave_resistance_odd_modes * wave_resistance_even_modes *
    (wave_resistance_odd_modes + wave_resistance_even_modes)**2) /
    ((wave_resistance_odd_modes + wave_resistance_even_modes * (number_segments - 1)) *
    (wave_resistance_even_modes + wave_resistance_odd_modes * (number_segments - 1)))))


def print_law() -> str:
    return print_expression(law)


@validate_input(wave_resistance_odd_modes_=wave_resistance_odd_modes,
    wave_resistance_even_modes_=wave_resistance_even_modes,
    number_segments_=number_segments)
@validate_output(wave_resistance)
def calculate_wave_resistance(wave_resistance_odd_modes_: Quantity,
    wave_resistance_even_modes_: Quantity, number_segments_: int) -> Quantity:
    result_expr = solve(law, wave_resistance, dict=True)[0][wave_resistance]
    result_expr = result_expr.subs({
        wave_resistance_odd_modes: wave_resistance_odd_modes_,
        wave_resistance_even_modes: wave_resistance_even_modes_,
        number_segments: number_segments_
    })
    return Quantity(result_expr)
