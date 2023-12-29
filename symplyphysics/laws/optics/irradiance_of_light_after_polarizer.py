from sympy import Eq, solve, cos
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
    angle_type,
)

# Description
## Malus's law states that the irradiance of linearly polarized light that passes through a polarizer
## is proportional to its initial intensity I0, transparency coefficient k of the polarizer
## and square cosine of the angle phi between the light's initial polarization direction and the axis
## of the polarizer

# Law: I = I0 * k * cos^2(phi)
## I is the irradiance of light passed through the polarizer
## I0 is the initial irradiance
## k is the polarizer's transparency coefficient
## phi is the angle between the initial polarization direction and the axis of the polarizer

irradiance_final = Symbol("irradiance_final", units.power / units.area)
irradiance_initial = Symbol("irradiance_initial", units.power / units.area)
transparency_coefficient = Symbol("transparency_coefficient", dimensionless)
polarization_angle = Symbol("polarization_angle", angle_type)

law = Eq(
    irradiance_final,
    irradiance_initial * transparency_coefficient * cos(polarization_angle) ** 2,
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    irradiance_initial_=irradiance_initial,
    transparency_coefficient_=transparency_coefficient,
    polarization_angle_=polarization_angle,
)
@validate_output(irradiance_final)
def calculate_irradiance(
    irradiance_initial_: Quantity,
    transparency_coefficient_: float,
    polarization_angle_: Quantity | float,
) -> Quantity:
    result_expr = solve(law, irradiance_final, dict=True)[0][irradiance_final]
    irradiance_applied = result_expr.subs({
            irradiance_initial: irradiance_initial_,
            transparency_coefficient: transparency_coefficient_,
            polarization_angle: polarization_angle_,
    })
    return Quantity(irradiance_applied)
