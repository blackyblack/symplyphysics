from sympy import Eq, solve, sqrt
from sympy.physics.units import magnetic_constant
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless, angle_type,)

## Description
## A coaxial waveguide is an electrical cable consisting of a central conductor and a shield arranged coaxially and separated
## by an insulating material or an air gap. It is used to transmit radio frequency electrical signals.
## The resistance formed by the surface of the metal sheath of the cable can be calculated by knowing the signal frequency,
## magnetic permeability and specific conductivity conductivity of the metal.

## Law is: Rs = sqrt(w * mu0 * mur / (2 * sig)), where
## Rs - surface resistance,
## w - angular frequency of signal,
## mu0 - magnetic constant,
## mur - relative permeability of the insulator material,
## sig - specific conductivity of conductor.

surface_resistance = Symbol("surface_resistance", units.impedance)

relative_permeability = Symbol("relative_permeability", dimensionless)
angular_frequency = Symbol("angular_frequency", angle_type / units.time)
specific_conductivity = Symbol("specific_conductivity", units.conductance / units.length)

law = Eq(surface_resistance, sqrt(angular_frequency * magnetic_constant * relative_permeability / (2 * specific_conductivity)))


def print_law() -> str:
    return print_expression(law)


@validate_input(relative_permeability_=relative_permeability, angular_frequency_=angular_frequency, specific_conductivity_=specific_conductivity)
@validate_output(surface_resistance)
def calculate_surface_resistance(relative_permeability_: float, angular_frequency_: Quantity,
    specific_conductivity_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, surface_resistance, dict=True)[0][surface_resistance]
    result_expr = result_velocity_expr.subs({
        relative_permeability: relative_permeability_,
        angular_frequency: angular_frequency_,
        specific_conductivity: specific_conductivity_
    })
    return Quantity(result_expr)
