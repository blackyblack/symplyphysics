from sympy import Eq, solve, sqrt
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output,)

## Description
## A rectangular waveguide is a rectangular metal waveguide capable of supporting waves propagating along it.
## There is a critical wavelength. Signals with a wavelength greater than the critical one are attenuated and
## do not propagate in the waveguide.
## The characteristic resistance of a wave is a value determined by the ratio of the transverse component
## of the electric field strength to the transverse component of the magnetic field strength of a traveling wave.

## Law is: Z = Z0 / sqrt(1 - (L / L1)^2), where
## Z - characteristic resistance of rectangular waveguide,
## Z0 - characteristic resistance of the material filling the waveguide,
## L - wavelength,
## L1 - critical wavelength.

resistance = Symbol("resistance", units.impedance)

resistance_of_medium = Symbol("resistance_of_medium", units.impedance)
wavelength = Symbol("wavelength", units.length)
critical_wavelength = Symbol("critical_wavelength", units.length)

law = Eq(resistance, resistance_of_medium / sqrt(1 - (wavelength / critical_wavelength)**2))


def print_law() -> str:
    return print_expression(law)


@validate_input(resistance_of_medium_=resistance_of_medium, wavelength_=wavelength, critical_wavelength_=critical_wavelength)
@validate_output(resistance)
def calculate_resistance(resistance_of_medium_: Quantity, wavelength_: Quantity,
    critical_wavelength_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, resistance, dict=True)[0][resistance]
    result_expr = result_velocity_expr.subs({
        resistance_of_medium: resistance_of_medium_,
        wavelength: wavelength_,
        critical_wavelength: critical_wavelength_
    })
    return Quantity(result_expr)
