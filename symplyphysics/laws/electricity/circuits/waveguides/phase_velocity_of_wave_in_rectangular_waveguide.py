from sympy import Eq, solve, sqrt
from sympy.physics.units import speed_of_light
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
## A rectangular waveguide is a rectangular metal waveguide capable of supporting waves propagating along it.
## There is a critical wavelength. Signals with a wavelength greater than the critical one are attenuated and
## do not propagate in the waveguide.

## Law is: Vph = c / sqrt(er * mur * (1 - (L / L1)^2)), where
## Vph - phase velocity of wave in rectangular waveguide,
## c - speed of light,
## er - relative permittivity of insulating material,
## mur - relative permeability of the insulating material,
## L - wavelength,
## L1 - critical wavelength.

phase_velocity = Symbol("phase_velocity", units.velocity)

relative_permittivity = Symbol("relative_permittivity", dimensionless)
relative_permeability = Symbol("relative_permeability", dimensionless)
wavelength = Symbol("wavelength", units.length)
critical_wavelength = Symbol("critical_wavelength", units.length)

law = Eq(
    phase_velocity, speed_of_light / sqrt(relative_permittivity * relative_permeability * (1 -
    (wavelength / critical_wavelength)**2)))


def print_law() -> str:
    return print_expression(law)


@validate_input(relative_permittivity_=relative_permittivity,
    relative_permeability_=relative_permeability,
    wavelength_=wavelength,
    critical_wavelength_=critical_wavelength)
@validate_output(phase_velocity)
def calculate_phase_velocity(relative_permittivity_: float, relative_permeability_: float,
    wavelength_: Quantity, critical_wavelength_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, phase_velocity, dict=True)[0][phase_velocity]
    result_expr = result_velocity_expr.subs({
        relative_permittivity: relative_permittivity_,
        relative_permeability: relative_permeability_,
        wavelength: wavelength_,
        critical_wavelength: critical_wavelength_
    })
    return Quantity(result_expr)
