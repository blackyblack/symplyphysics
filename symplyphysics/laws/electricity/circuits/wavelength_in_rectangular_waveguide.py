from sympy import Eq, solve, sqrt
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output,)

## Description
## A rectangular waveguide is a rectangular metal waveguide capable of supporting waves propagating along it.
## There is a critical wavelength. Signals with a wavelength greater than the critical one are attenuated and
## do not propagate in the waveguide.

## Law is: Lw = L / sqrt(1 - (L / L1)^2), where
## Lw - wavelength in a rectangular waveguide,
## L - wavelength,
## L1 - critical wavelength.

wavelength_waveguide = Symbol("wavelength_waveguide", units.length)

wavelength = Symbol("wavelength", units.length)
critical_wavelength = Symbol("critical_wavelength", units.length)

law = Eq(wavelength_waveguide, wavelength / sqrt(1 - (wavelength / critical_wavelength)**2))


def print_law() -> str:
    return print_expression(law)


@validate_input(wavelength_=wavelength, critical_wavelength_=critical_wavelength)
@validate_output(wavelength_waveguide)
def calculate_wavelength_waveguide(wavelength_: Quantity, critical_wavelength_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, wavelength_waveguide, dict=True)[0][wavelength_waveguide]
    result_expr = result_velocity_expr.subs({
        wavelength: wavelength_,
        critical_wavelength: critical_wavelength_
    })
    return Quantity(result_expr)
