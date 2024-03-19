from sympy import Eq, solve
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless)

## Description
## A coaxial waveguide is an electrical cable consisting of a central conductor and a shield arranged coaxially and separated
## by an insulating material or an air gap. It is used to transmit radio frequency electrical signals.
## The specific conductivity of a coaxial waveguide depends on the frequency of signal and the specific capacity of coaxial waveguide,
## as well as on the tangent of the dielectric loss angle of the insulator material.

## Law is: G = w * C * tan(d), where
## G - specific conductivity of coaxial waveguide,
## w - frequency of signal,
## C - specific capacity of coaxial waveguide,
## tan(d) - tangent of the dielectric loss angle of the insulator material.

specific_conductivity = Symbol("specific_conductivity", units.conductance / units.length)

frequency = Symbol("frequency", 1 / units.time)
specific_capacity = Symbol("specific_capacity", units.capacitance / units.length)
tangent_dielectric_loss_angle = Symbol("tangent_dielectric_loss_angle", dimensionless)

law = Eq(specific_conductivity, frequency * specific_capacity * tangent_dielectric_loss_angle)


def print_law() -> str:
    return print_expression(law)


@validate_input(frequency_=frequency, specific_capacity_=specific_capacity, tangent_dielectric_loss_angle_=tangent_dielectric_loss_angle)
@validate_output(specific_conductivity)
def calculate_specific_conductivity(frequency_: Quantity, specific_capacity_: Quantity,
    tangent_dielectric_loss_angle_: float) -> Quantity:
    result_velocity_expr = solve(law, specific_conductivity, dict=True)[0][specific_conductivity]
    result_expr = result_velocity_expr.subs({
        frequency: frequency_,
        specific_capacity: specific_capacity_,
        tangent_dielectric_loss_angle: tangent_dielectric_loss_angle_
    })
    return Quantity(result_expr)
