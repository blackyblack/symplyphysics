from sympy import Eq, solve, pi
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
## The attenuation coefficient shows how many times the transmitted signal weakens per unit length of the coaxial waveguide.

## Law is: ad = (pi / L) * tan(d) * Z / Z0, where
## ad - attenuation coefficient of waveguide,
## L - wavelength,
## tan(d) - tangent of the dielectric loss angle of the material filling the waveguide,
## Z - characteristic resistance of rectangular waveguide,
## Z0 - characteristic resistance of the material filling the waveguide.

attenuation_coefficient = Symbol("attenuation_coefficient", 1 / units.length)

resistance_of_waveguide = Symbol("resistance_of_waveguide", units.impedance)
resistance_of_medium = Symbol("resistance_of_medium", units.impedance)
wavelength = Symbol("wavelength", units.length)
tangent_dielectric_loss_angle = Symbol("tangent_dielectric_loss_angle", dimensionless)

law = Eq(attenuation_coefficient, (pi / wavelength) * tangent_dielectric_loss_angle *
    resistance_of_waveguide / resistance_of_medium)


def print_law() -> str:
    return print_expression(law)


@validate_input(resistance_of_waveguide_=resistance_of_waveguide,
    resistance_of_medium_=resistance_of_medium,
    wavelength_=wavelength,
    tangent_dielectric_loss_angle_=tangent_dielectric_loss_angle)
@validate_output(attenuation_coefficient)
def calculate_attenuation_coefficient(resistance_of_waveguide_: Quantity,
    resistance_of_medium_: Quantity, wavelength_: Quantity,
    tangent_dielectric_loss_angle_: float) -> Quantity:
    result_velocity_expr = solve(law, attenuation_coefficient,
        dict=True)[0][attenuation_coefficient]
    result_expr = result_velocity_expr.subs({
        resistance_of_waveguide: resistance_of_waveguide_,
        resistance_of_medium: resistance_of_medium_,
        wavelength: wavelength_,
        tangent_dielectric_loss_angle: tangent_dielectric_loss_angle_
    })
    return Quantity(result_expr)
