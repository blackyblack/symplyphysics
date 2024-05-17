from sympy import Eq, solve, sqrt
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
)

## Description
## The microstrip line is a dielectric substrate on which a metal strip is applied.
## When a wave propagates along a microstrip line, part of the field goes out, since the microstrip line does
## not have metal borders on all sides, unlike, for example, rectangular waveguides. Then imagine an environment
## in which the field will have the same magnitude as the field of a microstrip line. The permittivity of such a
## medium will be called the effective permittivity of the line.
## The attenuation coefficient shows how many times the transmitted signal weakens per unit length of the microstrip line.

## Law is: ad = 27.3 * (er / sqrt(ef)) * ((ef - 1) / (er - 1)) * tan(d) / L, where
## ad - attenuation coefficient of the microstrip line,
## er - relative permittivity of the dielectric substrate of the microstrip line,
## ef - effective permittivity of the microstrip line,
## tan(d) - tangent of the dielectric loss angle of the dielectric substrate of the microstrip line,
## L - wavelength in vacuum.

attenuation_coefficient = Symbol("attenuation_coefficient", 1 / units.length)

relative_permittivity = Symbol("relative_permittivity", dimensionless)
effective_permittivity = Symbol("effective_permittivity", dimensionless)
wavelength = Symbol("wavelength", units.length)
tangent_dielectric_loss_angle = Symbol("tangent_dielectric_loss_angle", dimensionless)

expression_1 = relative_permittivity / sqrt(effective_permittivity)
expression_2 = (effective_permittivity - 1) / (relative_permittivity - 1)
expression_3 = tangent_dielectric_loss_angle / wavelength
constant_coefficient = 27.3
law = Eq(attenuation_coefficient, constant_coefficient * expression_1 * expression_2 * expression_3)


@validate_input(relative_permittivity_=relative_permittivity,
    effective_permittivity_=effective_permittivity,
    wavelength_=wavelength,
    tangent_dielectric_loss_angle_=tangent_dielectric_loss_angle)
@validate_output(attenuation_coefficient)
def calculate_attenuation_coefficient(relative_permittivity_: float, effective_permittivity_: float,
    wavelength_: Quantity, tangent_dielectric_loss_angle_: float) -> Quantity:
    result_expr = solve(law, attenuation_coefficient, dict=True)[0][attenuation_coefficient]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        effective_permittivity: effective_permittivity_,
        wavelength: wavelength_,
        tangent_dielectric_loss_angle: tangent_dielectric_loss_angle_
    })
    return Quantity(result_expr)
