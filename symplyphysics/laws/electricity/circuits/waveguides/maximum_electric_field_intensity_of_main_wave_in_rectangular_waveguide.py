from sympy import Eq, solve, pi, sqrt
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless)

## Description
## A rectangular waveguide is a rectangular metal waveguide capable of supporting waves propagating along it.
## The main wave is a transverse electric wave with the first index equal to 1 and the second index equal to 0.
## The first index shows how many half-wave lengths fit horizontally across the cross section. The second index
## shows how many half-wave lengths fit vertically across the cross section.

## Law is: E = 2 * 120 * pi * a * H / (L * sqrt(er)), where
## E - maximum electric field intensity in waveguide,
## a - width of the waveguide cross section,
## H - magnetic field intensity in waveguide,
## L - wavelength,
## er - relative permittivity of the material filling the waveguide.

# Conditions:
# - the wave propagating in the waveguide must be the main wave;
# - the waveguide must be rectangular.

maximum_electric_intensity = Symbol("maximum_electric_intensity", units.voltage / units.length)

relative_permittivity = Symbol("relative_permittivity", dimensionless)
waveguide_width = Symbol("waveguide_width", units.length)
wavelength = Symbol("wavelength", units.length)
magnetic_intensity = Symbol("magnetic_intensity", units.current / units.length)

vacuum_impedance = Quantity(120 * units.impedance)

law = Eq(
    maximum_electric_intensity, 2 * vacuum_impedance * pi * waveguide_width * magnetic_intensity /
    (wavelength * sqrt(relative_permittivity)))


def print_law() -> str:
    return print_expression(law)


@validate_input(relative_permittivity_=relative_permittivity,
    waveguide_width_=waveguide_width,
    wavelength_=wavelength,
    magnetic_intensity_=magnetic_intensity)
@validate_output(maximum_electric_intensity)
def calculate_maximum_electric_intensity(relative_permittivity_: float, waveguide_width_: Quantity,
    wavelength_: Quantity, magnetic_intensity_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, maximum_electric_intensity,
        dict=True)[0][maximum_electric_intensity]
    result_expr = result_velocity_expr.subs({
        relative_permittivity: relative_permittivity_,
        waveguide_width: waveguide_width_,
        wavelength: wavelength_,
        magnetic_intensity: magnetic_intensity_
    })
    return Quantity(result_expr)
