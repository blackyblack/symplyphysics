from sympy import (
    Eq,
    solve,
    sqrt,
)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## A rectangular waveguide is a rectangular metal waveguide capable of supporting waves propagating along it.
## The main wave is a transverse electric wave with the first index equal to 1 and the second index equal to 0.
## The first index shows how many half-wave lengths fit horizontally across the cross section. The second index
## shows how many half-wave lengths fit vertically across the cross section.

## Law is: P = a * b * sqrt(1 - (L / (2 * a))^2) * E0^2 / (4 * Z0), where
## P - power carried by the waveguide,
## a - width of the waveguide cross section,
## b - height of the waveguide cross section,
## L - wavelength,
## E0 - maximum electric field intensity in the waveguide,
## Z0 - characteristic resistance of the material filling the waveguide.

# Conditions:
# - The wave propagating in the waveguide must be the main wave.

waveguide_power = Symbol("waveguide_power", units.power)

waveguide_width = Symbol("waveguide_width", units.length)
waveguide_height = Symbol("waveguide_height", units.length)
wavelength = Symbol("wavelength", units.length)
material_resistance = Symbol("material_resistance", units.impedance)
electric_intensity = Symbol("electric_intensity", units.voltage / units.length)

law = Eq(
    waveguide_power,
    waveguide_width * waveguide_height * sqrt(1 - (wavelength / (2 * waveguide_width))**2) *
    electric_intensity**2 / (4 * material_resistance))


def print_law() -> str:
    return print_expression(law)


@validate_input(waveguide_width_=waveguide_width,
    waveguide_height_=waveguide_height,
    wavelength_=wavelength,
    material_resistance_=material_resistance,
    electric_intensity_=electric_intensity)
@validate_output(waveguide_power)
def calculate_waveguide_power(waveguide_width_: Quantity, waveguide_height_: Quantity,
    wavelength_: Quantity, material_resistance_: Quantity,
    electric_intensity_: Quantity) -> Quantity:
    result_expr = solve(law, waveguide_power, dict=True)[0][waveguide_power]
    result_expr = result_expr.subs({
        waveguide_width: waveguide_width_,
        waveguide_height: waveguide_height_,
        wavelength: wavelength_,
        material_resistance: material_resistance_,
        electric_intensity: electric_intensity_
    })
    return Quantity(result_expr)
