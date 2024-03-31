from sympy import Eq, pi
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    angle_type,
)

# Description
## Wavenumber is the spatial frequency of a wave, measured in radians per unit distance.

# Law: k = 2*pi / lambda
## k - angular wavenumber
## lambda - wavelength

angular_wavenumber = Symbol("angular_wave_number", angle_type / units.length)
wavelength = Symbol("wavelength", units.length)

definition = Eq(angular_wavenumber, 2 * pi / wavelength)


def print_law() -> str:
    return print_expression(definition)


@validate_input(wavelength_=wavelength)
@validate_output(angular_wavenumber)
def calculate_wavenumber(wavelength_: Quantity) -> Quantity:
    result = definition.rhs.subs(wavelength, wavelength_)
    return Quantity(result)
