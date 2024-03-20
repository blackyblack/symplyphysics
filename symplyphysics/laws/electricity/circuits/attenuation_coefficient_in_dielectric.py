from sympy import Eq, solve, sqrt
from sympy.physics.units import electric_constant, magnetic_constant
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless, angle_type,)

## Description
## A coaxial waveguide is an electrical cable consisting of a central conductor and a shield arranged coaxially and separated
## by an insulating material or an air gap. It is used to transmit radio frequency electrical signals.
## The attenuation coefficient of a coaxial waveguide depends on the frequency of signal, as well as on the relative permittivity,
## the relative permeability and the dielectric loss angle of the insulator material.

## Law is: ad = (1 / 2) * w * sqrt(mu0 * mur * e0 * er) * tan(d), where
## ad - attenuation coefficient of coaxial waveguide,
## w - angular frequency of signal,
## e0 - electric constant,
## er - relative permittivity of insulating material,
## mu0 - magnetic constant,
## mur - relative permeability of the insulator material,
## tan(d) - tangent of the dielectric loss angle of the insulator material.

attenuation_coefficient = Symbol("attenuation_coefficient", 1 / units.length)

relative_permittivity = Symbol("relative_permittivity", dimensionless)
relative_permeability = Symbol("relative_permeability", dimensionless)
angular_frequency = Symbol("angular_frequency", angle_type / units.time)
tangent_dielectric_loss_angle = Symbol("tangent_dielectric_loss_angle", dimensionless)

law = Eq(attenuation_coefficient, (1 / 2) * angular_frequency * sqrt(magnetic_constant * relative_permeability * electric_constant * relative_permittivity) * tangent_dielectric_loss_angle)


def print_law() -> str:
    return print_expression(law)


@validate_input(relative_permittivity_=relative_permittivity, relative_permeability_=relative_permeability, angular_frequency_=angular_frequency, tangent_dielectric_loss_angle_=tangent_dielectric_loss_angle)
@validate_output(attenuation_coefficient)
def calculate_attenuation_coefficient(relative_permittivity_: float, relative_permeability_: float, angular_frequency_: Quantity,
    tangent_dielectric_loss_angle_: float) -> Quantity:
    result_velocity_expr = solve(law, attenuation_coefficient, dict=True)[0][attenuation_coefficient]
    result_expr = result_velocity_expr.subs({
        relative_permittivity: relative_permittivity_,
        relative_permeability: relative_permeability_,
        angular_frequency: angular_frequency_,
        tangent_dielectric_loss_angle: tangent_dielectric_loss_angle_
    })
    return Quantity(result_expr)
