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
## and odd modes are distributed. Knowing the coupling coefficient between the coupler segments, the wave resistance
## of the odd mode, as well as the number of coupler segments, it is possible to calculate the wave resistance for an even mode.
## https://habrastorage.org/r/w1560/getpro/habr/upload_files/054/d02/c8d/054d02c8d91c06425ae079d34b18ce15.jpeg

## Law is: Z = Zo * (C0 + sqrt(C0**2 + (1 - Co**2) * (N - 1)**2)) / ((N - 1) * (1 - C0)), where
## Z - wave resistance of the even mode,
## Zo - wave resistance of the odd mode,
## C0 - coupling factor between the coupler segments,
## N - the number of segments of the Lange coupler.

wave_resistance_even_modes = Symbol("wave_resistance_even_modes", units.impedance)

coupling_factor = Symbol("coupling_factor", dimensionless)
wave_resistance_odd_modes = Symbol("wave_resistance_odd_modes", units.impedance)
number_segments = Symbol("number_segments", dimensionless)

law = Eq(wave_resistance_even_modes, wave_resistance_odd_modes * (coupling_factor + sqrt(coupling_factor**2 + (1 - coupling_factor**2) * (number_segments - 1)**2)) / ((number_segments - 1) * (1 - coupling_factor)))

def print_law() -> str:
    return print_expression(law)


@validate_input(coupling_factor_=coupling_factor,
    wave_resistance_odd_modes_=wave_resistance_odd_modes,
    number_segments_=number_segments)
@validate_output(wave_resistance_even_modes)
def calculate_wave_resistance_even_modes(coupling_factor_: float, wave_resistance_odd_modes_: Quantity,
    number_segments_: int) -> Quantity:
    result_expr = solve(law, wave_resistance_even_modes, dict=True)[0][wave_resistance_even_modes]
    result_expr = result_expr.subs({
        coupling_factor: coupling_factor_,
        wave_resistance_odd_modes: wave_resistance_odd_modes_,
        number_segments: number_segments_
    })
    return Quantity(result_expr)
