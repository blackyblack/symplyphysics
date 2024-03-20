from sympy import Eq, solve
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless, angle_type,)

## Description
## A coaxial waveguide is an electrical cable consisting of a central conductor and a shield arranged coaxially and separated
## by an insulating material or an air gap. It is used to transmit radio frequency electrical signals.
## The specific conductivity of a coaxial waveguide depends on the frequency of signal and the specific capacitance of coaxial waveguide,
## as well as on the tangent of the dielectric loss angle of the insulator material.
## In general, the tangent of the dielectric loss angle is a characteristic of the material. The angle itself is usually not searched for or indicated.
## The tangent of this angle is found as the ratio of the active current to the reactive current. The value of this value for each material is in the tables
## in the public domain.

## Law is: G = w * C * tan(d), where
## G - specific conductivity of coaxial waveguide,
## w - angular frequency of signal,
## C - specific capacitance of coaxial waveguide,
## tan(d) - tangent of the dielectric loss angle of the insulator material.

specific_conductivity = Symbol("specific_conductivity", units.conductance / units.length)

angular_frequency = Symbol("angular_frequency", angle_type / units.time)
specific_capacitance = Symbol("specific_capacitance", units.capacitance / units.length)
tangent_dielectric_loss_angle = Symbol("tangent_dielectric_loss_angle", dimensionless)

law = Eq(specific_conductivity, angular_frequency * specific_capacitance * tangent_dielectric_loss_angle)


def print_law() -> str:
    return print_expression(law)


@validate_input(angular_frequency_=angular_frequency, specific_capacitance_=specific_capacitance, tangent_dielectric_loss_angle_=tangent_dielectric_loss_angle)
@validate_output(specific_conductivity)
def calculate_specific_conductivity(angular_frequency_: Quantity, specific_capacitance_: Quantity,
    tangent_dielectric_loss_angle_: float) -> Quantity:
    result_velocity_expr = solve(law, specific_conductivity, dict=True)[0][specific_conductivity]
    result_expr = result_velocity_expr.subs({
        angular_frequency: angular_frequency_,
        specific_capacitance: specific_capacitance_,
        tangent_dielectric_loss_angle: tangent_dielectric_loss_angle_
    })
    return Quantity(result_expr)
