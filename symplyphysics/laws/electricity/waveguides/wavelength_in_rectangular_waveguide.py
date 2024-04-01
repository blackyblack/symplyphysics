from sympy import Eq, solve, sqrt
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

## Description
## A rectangular waveguide is a rectangular metal waveguide capable of supporting waves propagating along it.
## There is a critical wavelength. Signals with a wavelength greater than the critical one are attenuated and
## do not propagate in the waveguide.

## Law is: Lw = L / sqrt(1 - (L / L1)^2), where
## Lw - wavelength in a rectangular waveguide,
## L - wavelength of a signal in a waveguide,
## L1 - critical wavelength.

waveguide_wavelength = Symbol("waveguide_wavelength", units.length)

signal_wavelength = Symbol("signal_wavelength", units.length)
critical_wavelength = Symbol("critical_wavelength", units.length)

law = Eq(waveguide_wavelength,
    signal_wavelength / sqrt(1 - (signal_wavelength / critical_wavelength)**2))


def print_law() -> str:
    return print_expression(law)


@validate_input(signal_wavelength_=signal_wavelength, critical_wavelength_=critical_wavelength)
@validate_output(waveguide_wavelength)
def calculate_waveguide_wavelength(signal_wavelength_: Quantity,
    critical_wavelength_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, waveguide_wavelength, dict=True)[0][waveguide_wavelength]
    result_expr = result_velocity_expr.subs({
        signal_wavelength: signal_wavelength_,
        critical_wavelength: critical_wavelength_
    })
    return Quantity(result_expr)
