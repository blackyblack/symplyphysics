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
    dimensionless,
)

# Description
## A coaxial waveguide is an electrical cable consisting of a central conductor and a shield arranged coaxially and separated
## by an insulating material or an air gap. It is used to transmit radio frequency electrical signals.
## The specific resistance of a coaxial waveguide depends on the radius of the outer conductor and the radius of the inner conductor,
## as well as on the relative permeability of the insulator material, frequency of signal and specific conductivity of conductor.
## The attenuation coefficient shows how many times the transmitted signal weakens per unit length of the coaxial waveguide.
## The first index shows how many half-wave lengths fit horizontally across the cross section. The second index
## shows how many half-wave lengths fit vertically across the cross section.

## Law is: am = (2 * Rs / (Z0 * b * sqrt(1 - (L / (2 * L1))^2))) * ((1 + b / a) * (L / (2 * L1))^2 + (1 - (L / (2 * L1))^2) * (b / a) * (b * n^2 / a + m^2) / ((b * n / a)^2 + m^2)), where
## am - attenuation coefficient in metal,
## Rs - surface resistance,
## m - first index,
## n - second index,
## a - height of the waveguide cross section,
## b - width of the waveguide cross section,
## Z0 - characteristic resistance of the material filling the waveguide,
## L - wavelength,
## L1 - critical wavelength.

# Conditions:
# - waves propagating in the waveguide must be transverse electric waves;
# - second index must be greater than or equal to 1.

attenuation_coefficient = Symbol("attenuation_coefficient", 1 / units.length)

surface_resistance = Symbol("surface_resistance", units.impedance)
first_index = Symbol("first_index", dimensionless)
second_index = Symbol("second_index", dimensionless)
width = Symbol("width", units.length)
height = Symbol("height", units.length)
resistance_of_medium = Symbol("resistance_of_medium", units.impedance)
signal_wavelength = Symbol("signal_wavelength", units.length)
critical_wavelength = Symbol("critical_wavelength", units.length)

expression_1 = 2 * surface_resistance / (resistance_of_medium * width * sqrt(1 -
    (signal_wavelength / (2 * critical_wavelength))**2))
expression_2 = (1 + width / height) * (signal_wavelength / (2 * critical_wavelength))**2
expression_3 = (width / height) * ((width / height) * second_index**2 + first_index**2)
expression_4 = (width * second_index / height)**2 + first_index**2

law = Eq(
    attenuation_coefficient,
    expression_1 * (expression_2 + (1 - (signal_wavelength /
    (2 * critical_wavelength))**2) * expression_3 / expression_4))


def print_law() -> str:
    return print_expression(law)


@validate_input(surface_resistance_=surface_resistance,
    first_index_=first_index,
    second_index_=second_index,
    width_=width,
    height_=height,
    resistance_of_medium_=resistance_of_medium,
    signal_wavelength_=signal_wavelength,
    critical_wavelength_=critical_wavelength)
@validate_output(attenuation_coefficient)
def calculate_attenuation_coefficient(surface_resistance_: Quantity, first_index_: float,
    second_index_: float, width_: Quantity, height_: Quantity, resistance_of_medium_: Quantity,
    signal_wavelength_: Quantity, critical_wavelength_: Quantity) -> Quantity:
    # pylint: disable=too-many-arguments, too-many-positional-arguments
    if second_index_ < 1:
        raise ValueError("The second index must be greater than or equal to 1")
    result_expr = solve(law, attenuation_coefficient, dict=True)[0][attenuation_coefficient]
    result_expr = result_expr.subs({
        surface_resistance: surface_resistance_,
        first_index: first_index_,
        second_index: second_index_,
        width: width_,
        height: height_,
        resistance_of_medium: resistance_of_medium_,
        signal_wavelength: signal_wavelength_,
        critical_wavelength: critical_wavelength_,
    })
    return Quantity(result_expr)
